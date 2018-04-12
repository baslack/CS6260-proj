#include <maya/MFnPlugin.h>
#include <maya/MTypeId.h> 
#include <maya/MMatrixArray.h>
#include <maya/MStringArray.h>
#include <maya/MPxSkinCluster.h> 
#include <maya/MItGeometry.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MDoubleArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MFnMatrixData.h>
#include <maya/MQuaternion.h>
#include <maya/MMatrix.h>
#include <maya/MTransformationMatrix.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MFnPointArrayData.h>

#include <maya/MPxGPUDeformer.h>
#include <maya/MPxGPUDeformer.h>
#include <maya/MEvaluationNode.h>
#include <maya/MGPUDeformerRegistry.h>
#include <maya/MOpenCLInfo.h>
#include <maya/MViewport2Renderer.h>
#include <clew/clew_cl.h>
#include <CL/cl.h>
#include <maya/MPxDeformerNode.h>

#include <vector>

//my additions
#define NAME "corSkinCluster"
#define DEFAULT_OMEGA 0.1
#include <maya/MProfiler.h>

class corSkinCluster : public MPxSkinCluster
{
public:
    static  void*   creator();
    
	static  MStatus initialize();

	virtual MStatus deform(MDataBlock    &block,
                           MItGeometry   &iter,
                           const MMatrix &mat,
                           unsigned int   multiIndex);
	
	virtual MStatus precomp(MDataBlock);
	
    static const MTypeId id;
	static const int _profileCategory;
	static MObject cor_valid;
	static MObject cor_ar;

private:
	static MStatus handle_to_doublearray(MArrayDataHandle&, MDoubleArray&);

	static MStatus similarity(MDoubleArray&, MDoubleArray&, int, double&);

	static MStatus qlerp(MQuaternion&, MQuaternion&, MQuaternion&);

	static const double omega;
};

class parallel_corSkinCluster : public corSkinCluster
{
public:
	virtual MStatus deform(MDataBlock    &block,
                           MItGeometry   &iter,
                           const MMatrix &mat,
                           unsigned int   multiIndex);
};

// const MTypeId parallel_corSkinCluster::id (0x0122573);

const MTypeId corSkinCluster::id( 0x22573 );

const double corSkinCluster::omega(DEFAULT_OMEGA);

const int corSkinCluster::_profileCategory(MProfiler::addCategory(NAME));

MObject corSkinCluster::cor_valid;

MObject corSkinCluster::cor_ar;

void* corSkinCluster::creator()
{
	void *node = new corSkinCluster();
	return node;
}

MStatus corSkinCluster::initialize()
{
	MGlobal::startErrorLogging("C:\\\\Users\\iam\\Desktop\\corSkinCluster_init_log");
	
	MStatus status = MStatus::kSuccess;
	
	MFnNumericAttribute numeric_fn;
	cor_valid = numeric_fn.create("Valid_Precomputation", "valid", MFnNumericData::kBoolean, 0.0,  &status);
	if (status != MS::kSuccess){
		MGlobal::doErrorLogEntry("corSkinCluster:  error setting up valid attr.\n");
		return status;
	}
	numeric_fn.setStorable(true);
	addAttribute(cor_valid);
	
	MPointArray temp_ar;
	MFnPointArrayData fn;
	MObject default_ar_obj = fn.create(temp_ar); 

	MFnTypedAttribute typed_fn;
	cor_ar = typed_fn.create("Centers_of_Rotation", "cor", MFnData::Type::kPointArray, default_ar_obj, &status);
	if (status != MS::kSuccess){
		MGlobal::doErrorLogEntry("corSkinCluster:  error setting up CoR point array attr.\n");
		return status;
	}
	numeric_fn.setStorable(true);
	addAttribute(cor_ar);

	MGlobal::closeErrorLog();

	return MStatus::kSuccess;
}

MStatus corSkinCluster::handle_to_doublearray(MArrayDataHandle &handle, MDoubleArray &vec){
	int count = 0;
	int max = handle.elementCount();
	double val = 0.0;
	for (count = 0; count < max; count++){
		if (handle.jumpToElement(count) != MStatus::kSuccess){
			vec[count] = 0.0;
		}else{
			vec[count] = handle.inputValue().asDouble();
		}
	}
	return MStatus::kSuccess;
}

MStatus corSkinCluster::similarity(MDoubleArray &weight_p, 
								   MDoubleArray &weight_v,
								   int number_of_transforms,
								   double &result)
{
	int j = 0;
	int k = 0;
	result = 0;
	double temp = 0;


	for (j = 0; j < number_of_transforms; j++){
		for (k = 0; k < number_of_transforms; k++){
			if (j != k){
				auto wpj = weight_p[j];
				auto wpk = weight_p[k];
				auto wvj = weight_v[j];
				auto wvk = weight_v[k];
				temp = wpj*wpk*wvj*wvk;
				temp *= exp(-(pow(wpj*wvk-wpk*wvj,2.0)/pow(omega,2.0)));
				result += temp;
			}
		}  // end k loop
	}  // end j loop

	return MStatus::kSuccess;
}

MStatus corSkinCluster::qlerp(MQuaternion& q_a, MQuaternion& q_b, MQuaternion& result){
	MStatus stat;
	double dot_product;
	double q_a_comp[4];
	double q_b_comp[4];
	stat = q_a.get(q_a_comp);
	if (stat != MStatus::kSuccess){
		std::cerr << "corSkinCluster::qlerp, unable to extract q_a" << std::endl;
		return MStatus::kFailure;
	}
	stat = q_b.get(q_b_comp);
	if (stat != MStatus::kSuccess){
		std::cerr << "corSkinCluster::qlerp, unable to extract q_b" << std::endl;
		return MStatus::kFailure;
	}
	dot_product = 0.0;
	for (int i = 0; i < 4; i++){
		dot_product += q_a_comp[i]*q_b_comp[i];
	}
	if (dot_product >= 0){
		result = q_a + q_b;
	}else{
		result = q_a - q_b;
	}
	return MStatus::kSuccess;
}


MStatus corSkinCluster::precomp(MDataBlock block)
{
	// MGlobal::startErrorLogging("C:\\\\Users\\iam\\Desktop\\corSkinCluster_precomp_log");
	MStatus stat;

	// load current cor_ar and clear it
	MDataHandle cor_arHandle = block.inputValue(cor_ar);
	MFnData::Type test = cor_arHandle.type();
	MObject cor_arData = cor_arHandle.data();
	if (cor_arData.hasFn(MFn::Type::kPointArrayData)){
	}else{
		return MS::kFailure;
	}
	MFnPointArrayData cor_arFn(cor_arData, &stat);
	if (stat.error()){
		stat.perror("corSkinCluster::precomp, unable to attached MFnPtAr to cor\n");
		return stat;
	}
	MPointArray cor_PA = cor_arFn.array();
	cor_PA.clear();

	// get mesh iterator
	MArrayDataHandle inputHandle = block.inputArrayValue(input);
	stat = inputHandle.jumpToArrayElement(0);
	if (stat != MS::kSuccess){
		return stat;
	}
	MDataHandle cor_IOGeoHandle = inputHandle.inputValue().child(inputGeom);
	MObject cor_IOGeoObj = cor_IOGeoHandle.asMesh();
	MItMeshPolygon T(cor_IOGeoObj, &stat);
	if (stat.error()){
		stat.perror("corSkinCluster::precomp, unable to get mesh iterator\n");
		return stat;
	}

	// get vertex iterator
	MItGeometry v_i(cor_IOGeoHandle, false, &stat);
	if (stat.error()){
		stat.perror("corSkinCluster::precomp, unable to get vertex iterator\n");
		return stat;
	}

	// weights
	MArrayDataHandle w_i = block.inputArrayValue(weightList);
	if ( w_i.elementCount() == 0 ) {
		// no weights - nothing to do
		return MStatus::kFailure;
	}

	// bone transforms
	MArrayDataHandle transformHandle = block.inputArrayValue(matrix);
	int num_transforms = transformHandle.elementCount();

	//calculate per triange area
	MDoubleArray tri_area;
	//calculate average vertex position
	MPointArray tri_avg_pos;
	//calculate average vertex weights
	std::vector<MDoubleArray> tri_avg_weights;
	int num_tris;
	int idx;
	MPoint alpha,beta,gamma;
	MPointArray tri_verts;
	MIntArray tri_idx;
	MVector beta_alpha;
	MVector gamma_alpha;
	MDoubleArray tri_avg_weight;
	
	// pre calc all the areas, average weights and positions
	while(!(T.isDone())){
		// each hit on the iterator returns a face
		// that face is made of multiple triangles
		// how many?
		stat = T.numTriangles(num_tris);
		if (stat == MStatus::kSuccess){
			// for each triangle
			for (idx = 0; idx < num_tris; idx++){
				// get the verts
				stat = T.getTriangle(idx, tri_verts, tri_idx, MSpace::kObject);  // switched this to world, from kObject
				if (stat.error()){
					stat.perror("corSkinCluster::precomp, unable to get triangle from iterator\n");
					return stat;
				}
				alpha = tri_verts[0];
				beta = tri_verts[1];
				gamma = tri_verts[2];
				// calc and store area of triangle
				beta_alpha = MVector(beta-alpha);
				gamma_alpha = MVector(gamma-alpha);
				stat = tri_area.append(((beta_alpha ^ gamma_alpha).length())*0.5);
				if (stat.error()){
					stat.perror("corskinCluster::precomp, unable to append area\n");
					return stat;
				}
				
				// calc and store average vertex position
				stat = tri_avg_pos.append((alpha+beta+gamma)/3);
				if (stat.error()){
					stat.perror("corskinCluster::precomp, unable to apped average position\n");
					return stat;
				}
				
				// calc and store avg weights

				// get alpha weights
				stat = w_i.jumpToElement(tri_idx[0]);
				if (stat.error()){
					stat.perror("corSkinCluster::precomp, unable to get weights for alpha.\n");
					return stat;
				}
				MArrayDataHandle alpha_weightsHandle = w_i.inputValue().child(weights);

				// get beta weights
				stat = w_i.jumpToElement(tri_idx[1]);
				if (stat.error()){
					stat.perror("corSkinCluster::precomp, unable to get weights for beta.\n");
					return stat;
				}
				MArrayDataHandle beta_weightsHandle = w_i.inputValue().child(weights);

				// get gamma weights
				stat = w_i.jumpToElement(tri_idx[2]);
				if (stat.error()){
					stat.perror("corSkinCluster::precomp, unable to get weights for gamma.\n");
					return stat;
				}
				MArrayDataHandle gamma_weightsHandle = w_i.inputValue().child(weights);

				double a, b, c;
				tri_avg_weight.clear();
				for (int i = 0; i < num_transforms; i++){
					// average and store weights
					// get ith weight for alpha
					stat = alpha_weightsHandle.jumpToElement(i);
					if (stat == MStatus::kSuccess){
						a = alpha_weightsHandle.inputValue().asDouble();
					}else{
						a = 0.0;
					}
					// get ith weight for beta
					stat = beta_weightsHandle.jumpToElement(i);
					if (stat == MStatus::kSuccess){
						b = beta_weightsHandle.inputValue().asDouble();
					}else{
						b = 0.0;
					}
					// get ith weight for gamma
					stat = gamma_weightsHandle.jumpToElement(i);
					if (stat == MStatus::kSuccess){
						c = gamma_weightsHandle.inputValue().asDouble();
					}else{
						c = 0.0;
					}
					stat = tri_avg_weight.append((a + b + c)/3.0);
					if (stat.error()){
						stat.perror("corSkinCluster::precomp, unable to add average weight to array.\n");
						return stat;
					}
				} // end for
				tri_avg_weights.push_back(tri_avg_weight);
			} // end for
		}else{ // if num triangles fail
			stat.perror("corSkinCluster::precomp, face has no triangles?\n");
			return stat;
		} //end else
		T.next();  // next face
	} // end while, weights averaged, vertex positions averaged

	w_i.jumpToElement(0);
	v_i.reset();

	MPoint cor;
	double s, lower;
	MPoint upper;
	num_tris = tri_area.length();
	MDoubleArray vertex_weights;

	// for each point in the iterator
	while(!(v_i.isDone())){
		// calculate the COR

		// get the vertex weights in a double array
		vertex_weights.clear();
		MArrayDataHandle vertex_weights_handle = w_i.inputValue().child(weights);
		for (int i = 0; i < num_transforms; i++){
			stat = vertex_weights_handle.jumpToElement(i);
			if (stat.error()){
				vertex_weights.append(0.0);
			}else{
				vertex_weights.append(vertex_weights_handle.inputValue().asDouble());
			}
		}

		s = 0.0;
		upper = MPoint(0.0,0.0,0.0);
		lower = 0.0;
		// for each triangle
		for (int i = 0; i < num_tris; i++){
			stat = similarity(vertex_weights, tri_avg_weights[i], num_transforms, s);
			upper += tri_avg_pos[i]*s*tri_area[i];
			lower += s*tri_area[i];
		}
		if (lower > 0){
			cor = upper/lower;
		}else{
			cor = MPoint(0.0,0.0,0.0);
		}
		cor_PA.append(cor);
		
		// iterate the loop
		v_i.next();
		w_i.next();
	} // end while

	// put the computed point array back on the attribute
	cor_arFn.set(cor_PA);

	MGlobal::closeErrorLog();

	return MStatus::kSuccess;
}

MStatus corSkinCluster::deform( MDataBlock& block,
                      MItGeometry& iter,
                      const MMatrix& /*m*/,
                      unsigned int multiIndex)
//
// Method: deform
//
// Description:   Deforms the point with a simple smooth skinning algorithm
//
// Arguments:
//   block      : the datablock of the node
//   iter       : an iterator for the geometry to be deformed
//   m          : matrix to transform the point into world space
//   multiIndex : the index of the geometry that we are deforming
//
//
{
	MGlobal::startErrorLogging("C:\\\\Users\\iam\\Desktop\\corSkinCluster_deform_log");

	MStatus returnStatus;
    
	// get the influence transforms
	//
	MArrayDataHandle transformsHandle = block.inputArrayValue( matrix );
	int numTransforms = transformsHandle.elementCount();
	if ( numTransforms == 0 ) {
		return MS::kSuccess;
	}

	int precomp_event = MProfiler::eventBegin(corSkinCluster::_profileCategory, MProfiler::kColorG_L1, "corSkinCluster: precomp");

	// insert precomp test here
	MDataHandle validHandle = block.inputValue(cor_valid, &returnStatus);
	if (returnStatus != MS::kSuccess){
		MGlobal::doErrorLogEntry("corSkinCluster::deform, unable to get valid handle off datablock\n");
		return returnStatus;
	}
	if (!validHandle.asBool()){
		returnStatus = precomp(block);
		if (returnStatus != MS::kSuccess){
			MGlobal::doErrorLogEntry("corSkinCluster::deform, precomp returned error");
			return returnStatus;
		}
		validHandle.setBool(true);
	}

	// get the CORs
	MDataHandle cor_arHandle = block.inputValue(cor_ar, &returnStatus);
	if (returnStatus != MS::kSuccess){
		MGlobal::doErrorLogEntry("corSkinCluster::deform, unable to get cor_ar handle off datablock\n");
		return returnStatus;
	}
	MObject cor_arData = cor_arHandle.data();
	MFnPointArrayData cor_arFn(cor_arData, &returnStatus);
	if (returnStatus != MS::kSuccess){
		MGlobal::doErrorLogEntry("corSkinCluster::deform, unable to attach function set for cor ar attr obj\n");
		return returnStatus;
	}
	MPointArray cor_PA = cor_arFn.array();

	MProfiler::eventEnd(precomp_event);

	MMatrixArray transforms;
	for ( int i=0; i<numTransforms; ++i ) {
		transforms.append( MFnMatrixData( transformsHandle.inputValue().data() ).matrix() );
		transformsHandle.next();
	}

	MArrayDataHandle bindHandle = block.inputArrayValue( bindPreMatrix );
	if ( bindHandle.elementCount() > 0 ) {
		for ( int i=0; i<numTransforms; ++i ) {
			transforms[i] = MFnMatrixData( bindHandle.inputValue().data() ).matrix() * transforms[i];
			bindHandle.next();
		}
	}

	// get the unit quaternions for the rotations of the matricies
	std::vector<MQuaternion> q_j;
	for (int j=0; j<numTransforms; ++j ) {
		MTransformationMatrix temp_tm(transforms[j]);
		MQuaternion q = temp_tm.rotation().normalizeIt();
		q_j.push_back(q);
	}

	MArrayDataHandle weightListHandle = block.inputArrayValue( weightList );
	if ( weightListHandle.elementCount() == 0 ) {
		// no weights - nothing to do
		return MS::kSuccess;
	}

    // Iterate through each point in the geometry.
    //
    
	int skin_event = MProfiler::eventBegin(corSkinCluster::_profileCategory, MProfiler::kColorG_L1, "corSkinCluster: deform loop");
	
	for (int i = 0; !iter.isDone(); i++){
		MPoint v_i = iter.position();
		MPoint vprime_i;

		MArrayDataHandle weightsHandle = weightListHandle.inputValue().child(weights);

		MQuaternion q(0.0,0.0,0.0,0.0);
		MQuaternion result;
		double w_j;

		// find weighted q for v_i
		for (int j = 0; j < numTransforms;  j++){
			MStatus stat;
			stat = weightsHandle.jumpToElement(j);
			if (stat == MStatus::kSuccess){
				w_j = weightsHandle.inputValue().asDouble();
			}else{
				w_j = 0.0;
			}
			double q_comp[4];
			q_j[j].get(q_comp);
			stat = qlerp(q,w_j*q_j[j], q);
		}

		q = q.normalizeIt();
		MMatrix R = q.asMatrix();

		//LBS Matrix
		MMatrix R_t_prime = MMatrix::identity;  // init to identity
		for (int j = 0; j < numTransforms; j++){
			MTransformationMatrix tm(transforms[j]); 
			MTransformationMatrix temp = tm.asRotateMatrix();
			temp.setTranslation(tm.getTranslation(MSpace::kWorld, NULL), MSpace::kWorld);  // switched from kWorld
			
			MStatus stat;
			stat = weightsHandle.jumpToElement(j);
			if (stat == MStatus::kSuccess){
				w_j = weightsHandle.inputValue().asDouble();
			}else{
				w_j = 0.0;
			}
			if (j == 0){
				R_t_prime = w_j*temp.asMatrix();
			}else{
				R_t_prime += (w_j * temp.asMatrix());
			}
		}

		MTransformationMatrix R_t_prime_tm(R_t_prime);

		MVector t = (cor_PA[i] * R_t_prime_tm.asRotateMatrix()) + R_t_prime_tm.getTranslation(MSpace::kWorld, NULL) - (cor_PA[i] * R);  // switch from kWorld
		vprime_i = v_i * R + t;
		iter.setPosition(vprime_i);
		weightListHandle.next();
		iter.next();
	} //end skinning loop

	/**
	for ( ; !iter.isDone(); iter.next()) {
        MPoint pt = iter.position();
		MPoint skinned;

		// get the weights for this point
		MArrayDataHandle weightsHandle = weightListHandle.inputValue().child( weights );

		// compute the skinning
		for ( int i=0; i<numTransforms; ++i ) {
			if ( MS::kSuccess == weightsHandle.jumpToElement( i ) ) {
				skinned += ( pt * transforms[i] ) * weightsHandle.inputValue().asDouble();
			}
		}
        
		// Set the final position.
		iter.setPosition( skinned );

		// advance the weight list handle
		weightListHandle.next();
    }
	**/

	MProfiler::eventEnd(skin_event);
	MGlobal::closeErrorLog();
    return returnStatus;
}



// GPU Override

class corSkinGPUDeformer : public MPxGPUDeformer
{
public:
    // Virtual methods from MPxGPUDeformer
    corSkinGPUDeformer();
    virtual ~corSkinGPUDeformer();
    
    virtual MPxGPUDeformer::DeformerStatus evaluate(MDataBlock& block, const MEvaluationNode&, const MPlug& plug, unsigned int numElements, const MAutoCLMem, const MAutoCLEvent, MAutoCLMem, MAutoCLEvent&);
    virtual void terminate();
    static MGPUDeformerRegistrationInfo* getGPUDeformerInfo();
    static bool validateNodeInGraph(MDataBlock& block, const MEvaluationNode&, const MPlug& plug, MStringArray* messages);
    static bool validateNodeValues(MDataBlock& block, const MEvaluationNode&, const MPlug& plug, MStringArray* messages);
private:
    // helper methods
    void extractWeightArray(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug);
    void extractCoR(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug);
	void extractTMnQ(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug);
    // Storage for data on the GPU
    MAutoCLMem fCLWeights;
    MAutoCLMem fCoR;
	MAutoCLMem fTM;
	MAutoCLMem fQ;
	unsigned int fNumTransforms;
    unsigned int fNumElements;
    // Kernel
    MAutoCLKernel fKernel;
};

class corSkinNodeGPUDeformerInfo : public MGPUDeformerRegistrationInfo
{
public:
    corSkinNodeGPUDeformerInfo(){}
    virtual ~corSkinNodeGPUDeformerInfo(){}
    virtual MPxGPUDeformer* createGPUDeformer()
    {
        return new corSkinGPUDeformer();
    }
    
    virtual bool validateNodeInGraph(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages)
    {
        return corSkinGPUDeformer::validateNodeInGraph(block, evaluationNode, plug, messages);
    }
    virtual bool validateNodeValues(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages)
    {
        return corSkinGPUDeformer::validateNodeValues(block, evaluationNode, plug, messages);
    }
};

MGPUDeformerRegistrationInfo* corSkinGPUDeformer::getGPUDeformerInfo()
{
    static corSkinNodeGPUDeformerInfo theOne;
    return &theOne;
}

corSkinGPUDeformer::corSkinGPUDeformer()
{
}

corSkinGPUDeformer::~corSkinGPUDeformer()
{
    terminate();
}

// static
bool corSkinGPUDeformer::validateNodeInGraph(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages)
{
    return true;
}

// static
bool corSkinGPUDeformer::validateNodeValues(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages)
{   
    MObject node = plug.node();
    MFnDependencyNode fnNode(node);
 	MPlug validPlug(node, corSkinCluster::cor_valid);
    MDataHandle validData;
    validPlug.getValue(validData);
	if (!validData.asBool())
    {
        MOpenCLInfo::appendMessage(messages, "corSKin %s not supported by deformer evaluator because precomputation is invalid.", fnNode.name().asChar());
        return false;
    }
    return true;
}

MPxGPUDeformer::DeformerStatus corSkinGPUDeformer::evaluate(
    MDataBlock& block,                          // data block for "this" node
    const MEvaluationNode& evaluationNode,      // evaluation node representing "this" node
    const MPlug& plug,                          // the multi index we're working on.  There will be a separate instance created per multi index
    unsigned int numElements,                   // the number of float3 elements in inputBuffer and outputBuffer
    const MAutoCLMem inputBuffer,               // the input positions we are going to deform
    const MAutoCLEvent inputEvent,              // the input event we need to wait for before we start reading the input positions
    MAutoCLMem outputBuffer,                    // the output positions we should write to.  This may or may not be the same buffer as inputBuffer.
    MAutoCLEvent& outputEvent)                  // the event a downstream deformer will wait for before reading from output buffer
{
    // evaluate has two main pieces of work.  I need to transfer any data I care about onto the GPU, and I need to run my OpenCL Kernel.
    // First, transfer the data.  offset has two pieces of data I need to transfer to the GPU, the weight array and the offset matrix.
    // I don't need to transfer down the input position buffer, that is already handled by the deformer evaluator, the points are in inputBuffer.
    fNumElements = numElements;
    MObject node = plug.node();
    extractWeightArray(block, evaluationNode, plug);
    extractCoR(block, evaluationNode, plug);
	extractTMnQ(block, evaluationNode, plug);
    // Now that all the data we care about is on the GPU, setup and run the OpenCL Kernel
    if (!fKernel.get())
    {
        // char *maya_location = getenv( "MAYA_LOCATION" );
        // MString openCLKernelFile(maya_location);
        // openCLKernelFile +="/../Extra/devkitBase/devkit/plug-ins/offsetNode/offset.cl";
        MString openCLKernelFile("C:\\Users\\iam\\Documents\\maya\\2017\\plug-ins\\corSkin.cl");
		MString openCLKernelName("corSkin");
        fKernel = MOpenCLInfo::getOpenCLKernel(openCLKernelFile, openCLKernelName);
        if (!fKernel) return MPxGPUDeformer::kDeformerFailure;
    }
    cl_int err = CL_SUCCESS;
    
    // Set all of our kernel parameters.  Input buffer and output buffer may be changing every frame
    // so always set them.
    unsigned int parameterId = 0;
    err = clSetKernelArg(fKernel.get(), parameterId++, sizeof(cl_mem), (void*)outputBuffer.getReadOnlyRef());
    MOpenCLInfo::checkCLErrorStatus(err);
    err = clSetKernelArg(fKernel.get(), parameterId++, sizeof(cl_mem), (void*)inputBuffer.getReadOnlyRef());
    MOpenCLInfo::checkCLErrorStatus(err);
    err = clSetKernelArg(fKernel.get(), parameterId++, sizeof(cl_mem), (void*)fCLWeights.getReadOnlyRef());
    MOpenCLInfo::checkCLErrorStatus(err);
    err = clSetKernelArg(fKernel.get(), parameterId++, sizeof(cl_mem), (void*)fCoR.getReadOnlyRef());
    MOpenCLInfo::checkCLErrorStatus(err);
	err = clSetKernelArg(fKernel.get(), parameterId++, sizeof(cl_mem), (void*)fTM.getReadOnlyRef());
    MOpenCLInfo::checkCLErrorStatus(err);
    err = clSetKernelArg(fKernel.get(), parameterId++, sizeof(cl_mem), (void*)fQ.getReadOnlyRef());
    MOpenCLInfo::checkCLErrorStatus(err);
    err = clSetKernelArg(fKernel.get(), parameterId++, sizeof(cl_uint), (void*)&fNumElements);
    MOpenCLInfo::checkCLErrorStatus(err);

    // Figure out a good work group size for our kernel.
    size_t workGroupSize;
    size_t retSize;
    err = clGetKernelWorkGroupInfo(
        fKernel.get(),
        MOpenCLInfo::getOpenCLDeviceId(),
        CL_KERNEL_WORK_GROUP_SIZE,
        sizeof(size_t),
        &workGroupSize,
        &retSize);
    MOpenCLInfo::checkCLErrorStatus(err);
    size_t localWorkSize = 256;
    if (retSize > 0) localWorkSize = workGroupSize;
    size_t globalWorkSize = (localWorkSize - fNumElements % localWorkSize) + fNumElements; // global work size must be a multiple of localWorkSize
    // set up our input events.  The input event could be NULL, in that case we need to pass
    // slightly different parameters into clEnqueueNDRangeKernel
    unsigned int numInputEvents = 0;
    if (inputEvent.get())
    {
        numInputEvents = 1;
    }
    // run the kernel
    err = clEnqueueNDRangeKernel(
        MOpenCLInfo::getMayaDefaultOpenCLCommandQueue(),
        fKernel.get(),
        1,
        NULL,
        &globalWorkSize,
        &localWorkSize,
        numInputEvents,
        numInputEvents ? inputEvent.getReadOnlyRef() : 0,
        outputEvent.getReferenceForAssignment() );
    MOpenCLInfo::checkCLErrorStatus(err);
    return MPxGPUDeformer::kDeformerSuccess;
}

void corSkinGPUDeformer::terminate()
{
    MHWRender::MRenderer::theRenderer()->releaseGPUMemory(fNumElements*sizeof(float));
    fCLWeights.reset();
    fCoR.reset();
	fTM.reset();
	fQ.reset();
    MOpenCLInfo::releaseOpenCLKernel(fKernel);
    fKernel.reset();
}

void corSkinGPUDeformer::extractWeightArray(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug)
{
    // if we've already got a weight array and it is not changing then don't bother copying it
    // to the GPU again
    MStatus status;
    // Note that right now dirtyPlugExists takes an attribute, so if any element in the multi is changing we think it is dirty...
    // To avoid false dirty issues here you'd need to only use one element of the MPxDeformerNode::input multi attribute for each
    // offset node.
	if ((fCLWeights.get() && !evaluationNode.dirtyPlugExists(corSkinCluster::weightList, &status)) || !status)
    {
        return;
    }

	// what do we need to do
	// get the weight list
	// for each element of the weight list, push each of the weights

	std::vector<float> temp;

	MArrayDataHandle transformsHandle = block.inputArrayValue( corSkinCluster::matrix );
	int numTransforms = transformsHandle.elementCount();
	if ( numTransforms == 0 ) return;

	MStatus stat;
	MArrayDataHandle weightListHandle = block.inputArrayValue(corSkinCluster::weightList, &stat);
	if (stat.error()) return;
	for (unsigned int i = 0; i < fNumElements; i++){
		MArrayDataHandle weightsHandle = weightListHandle.inputValue().child(corSkinCluster::weights);
		for(int j = 0; j < numTransforms; j ++){
			stat = weightsHandle.jumpToElement(j);
			if (stat.error()){
				temp.push_back(0.0);
			}else{
				temp.push_back(weightsHandle.inputValue().asFloat());
			}
		}
		weightListHandle.next();
	}

	// weights all in temp
	
    // Two possibilities, we could be updating an existing OpenCL buffer or allocating a new one.
    cl_int err = CL_SUCCESS;
    if (!fCLWeights.get())
    {
        MHWRender::MRenderer::theRenderer()->holdGPUMemory(fNumElements*numTransforms*sizeof(float));
        fCLWeights.attach(clCreateBuffer(MOpenCLInfo::getOpenCLContext(), CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, fNumElements * numTransforms * sizeof(float), (void*)&temp[0], &err));
    }
    else
    {
        // I use a blocking write here, non-blocking could be faster...  need to manage the lifetime of temp, and have the kernel wait until the write finishes before running
        // I'm also assuming that the weight buffer is not growing.
        err = clEnqueueWriteBuffer(MOpenCLInfo::getMayaDefaultOpenCLCommandQueue(), fCLWeights.get(), CL_TRUE, 0, fNumElements * numTransforms * sizeof(float), (void*)&temp[0], 0, NULL, NULL);
    }
}

void corSkinGPUDeformer::extractCoR(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug)
{
	// get CoRs
	MStatus stat;

	// check for existing CoR in buffer
	/*
	if ((fCoR.get() && !evaluationNode.dirtyPlugExists(corSkinCluster::cor_ar, &stat)) || !stat)
    {
        return;
    }
	*/

	MDataHandle cor_arHandle = block.inputValue(corSkinCluster::cor_ar, &stat);
	if (stat!= MS::kSuccess){
		return;
	}
	MObject cor_arData = cor_arHandle.data();
	MFnPointArrayData cor_arFn(cor_arData, &stat);
	if (stat!= MS::kSuccess){
		return;
	}
	MPointArray cor_PA = cor_arFn.array();

	std::vector<float>temp;

	for (unsigned int i = 0; i < cor_PA.length(); i++){
		temp.push_back((float)cor_PA[i].x);
		temp.push_back((float)cor_PA[i].y);
		temp.push_back((float)cor_PA[i].z);
		temp.push_back((float)cor_PA[i].w);
	}

	// all CoR in temp

	cl_int err = CL_SUCCESS;
    if (!fCLWeights.get())
    {
        MHWRender::MRenderer::theRenderer()->holdGPUMemory(fNumElements * 4 * sizeof(float));
        fCLWeights.attach(clCreateBuffer(MOpenCLInfo::getOpenCLContext(), CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, fNumElements * 4 * sizeof(float), (void*)&temp[0], &err));
    }
    else
    {
        // I use a blocking write here, non-blocking could be faster...  need to manage the lifetime of temp, and have the kernel wait until the write finishes before running
        // I'm also assuming that the weight buffer is not growing.
        err = clEnqueueWriteBuffer(MOpenCLInfo::getMayaDefaultOpenCLCommandQueue(), fCLWeights.get(), CL_TRUE, 0, fNumElements * 4 * sizeof(float), (void*)&temp[0], 0, NULL, NULL);
    }
}

void corSkinGPUDeformer::extractTMnQ(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug)
{	
	// I pass the offset matrix to OpenCL using a buffer as well.  I also send down the inverse matrix to avoid calculating it many times on the GPU
	MStatus status;
	if ((fTM.get() && !evaluationNode.dirtyPlugExists(corSkinCluster::matrix , &status)) || !status)
	{
		return;
	}

	MArrayDataHandle transformsHandle = block.inputArrayValue( corSkinCluster::matrix );
	int numTransforms = transformsHandle.elementCount();
	if ( numTransforms == 0 ) {
		return;
	}
	
	MMatrixArray transforms;
	for ( int i=0; i<numTransforms; ++i ) {
		transforms.append( MFnMatrixData( transformsHandle.inputValue().data() ).matrix() );
		transformsHandle.next();
	}

	MArrayDataHandle bindHandle = block.inputArrayValue( corSkinCluster::bindPreMatrix );
	if ( bindHandle.elementCount() > 0 ) {
		for ( int i=0; i<numTransforms; ++i ) {
			transforms[i] = MFnMatrixData( bindHandle.inputValue().data() ).matrix() * transforms[i];
			bindHandle.next();
		}
	}

	// openCL needs matrices in transpose orientation
	// we also need to split off the scaling and shearing
	// as that the algorithm doesn't use that as of yet

	MMatrixArray transforms_transpose;
	MMatrix temp;
	MVector temp_trans;
	MTransformationMatrix temp_tm;
	for (int i = 0; i < numTransforms; i++){
		temp_tm = MTransformationMatrix(transforms[i]);
		temp_trans = temp_tm.getTranslation(MSpace::kWorld);
		temp = temp_tm.asRotateMatrix();
		temp_tm = MTransformationMatrix(temp);
		temp_tm.setTranslation(temp_trans, MSpace::kWorld);
		temp = temp_tm.asMatrix();
		temp = temp.transpose();
		transforms_transpose.append(temp);
	}

	// now we need to pack those transforms up into something
	// openCL can understand, i.e. and array of floats

	std::vector<float> conv_transforms;
	for ( int i = 0; i < numTransforms; i++){
		for ( int r = 0; r < 4; r++){
			for ( int c = 0; c < 4; c++){
				conv_transforms.push_back((float)transforms_transpose[i](r,c));
			}
		}
	}

	// now we prep that group for openCL
	cl_int err = CL_SUCCESS;
	if (!fTM.get())
	{
		fTM.attach(clCreateBuffer(MOpenCLInfo::getOpenCLContext(), CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, conv_transforms.size() * sizeof(float), (void*)&conv_transforms[0], &err));
	}
	else
	{
		// I use a blocking write here, non-blocking could be faster...  need to manage the lifetime of temp, and have the kernel wait until the write finishes before running
		err = clEnqueueWriteBuffer(MOpenCLInfo::getMayaDefaultOpenCLCommandQueue(), fTM.get(), CL_TRUE, 0, conv_transforms.size() * sizeof(float), (void*)&conv_transforms[0], 0, NULL, NULL);
	}

	// now, as that we already have the transforms here, 
	// let's get the quaternions while we're at it

	// get the unit quaternions for the rotations of the matricies
	std::vector<MQuaternion> q_j;
	for (int j=0; j<numTransforms; ++j ) {
		MTransformationMatrix temp_tm(transforms[j]);
		MQuaternion q = temp_tm.rotation().normalizeIt();
		q_j.push_back(q);
	}

	// now that we have the unit quaternions we need
	// to back them up in a fashion that openCL can handle them

	std::vector<float>conv_q;

	for (int j = 0; j < numTransforms;  ++j){
		conv_q.push_back((float)q_j[j].x);
		conv_q.push_back((float)q_j[j].y);
		conv_q.push_back((float)q_j[j].z);
		conv_q.push_back((float)q_j[j].w);
	}

	// pipe them to the device
	cl_int err = CL_SUCCESS;
	if (!fQ.get())
	{
		fQ.attach(clCreateBuffer(MOpenCLInfo::getOpenCLContext(), CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, conv_q.size() * sizeof(float), (void*)&conv_q[0], &err));
	}
	else
	{
		// I use a blocking write here, non-blocking could be faster...  need to manage the lifetime of temp, and have the kernel wait until the write finishes before running
		err = clEnqueueWriteBuffer(MOpenCLInfo::getMayaDefaultOpenCLCommandQueue(), fQ.get(), CL_TRUE, 0, conv_q.size() * sizeof(float), (void*)&conv_q[0], 0, NULL, NULL);
	}

	fNumTransforms = (unsigned int)numTransforms;
}


// standard initialization procedures
//

MStatus initializePlugin( MObject obj )
{
    MStatus result;

    MFnPlugin plugin( obj, "Benjamin Slack", "0.0", "Any");
    result = plugin.registerNode(
        "corSkinCluster" ,
        corSkinCluster::id ,
        &corSkinCluster::creator ,
        &corSkinCluster::initialize ,
        MPxNode::kSkinCluster
        );

    return result;
}

MStatus uninitializePlugin( MObject obj )
{
    MStatus result;

    MFnPlugin plugin( obj );
    result = plugin.deregisterNode( corSkinCluster::id );

    return result;
}
