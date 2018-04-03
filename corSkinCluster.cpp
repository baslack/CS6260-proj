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
	// cor_ar = typed_fn.create("Centers_of_Rotation", "cor", MFnData::Type::kPointArray, MObject::kNullObj, &status);
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
	double cross_product;
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
	cross_product = MVector(MPoint(q_a_comp))*MVector(MPoint(q_b_comp));
	if (cross_product >= 0){
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
		MQuaternion q = temp_tm.rotation();
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
