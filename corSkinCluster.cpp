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

    // Deformation function
    //
    virtual MStatus deform(MDataBlock    &block,
                           MItGeometry   &iter,
                           const MMatrix &mat,
                           unsigned int   multiIndex);
	
	virtual MStatus precomp(MDataBlock);
	
    static const MTypeId id;
	static const int _profileCategory;
private:
	MStatus handle_to_doublearray(MArrayDataHandle&, MDoubleArray&);

	virtual MStatus similarity(MDoubleArray&, MDoubleArray&, int, double&);

	static MStatus qlerp(MQuaternion&, MQuaternion&, MQuaternion&);

	static const double omega;

	virtual bool get_cor_valid();

	virtual void set_cor_valid(bool);
	
	bool cor_valid;

	//std::vector<double3> cor_vec;

	MPointArray cor_ar;

};

const MTypeId corSkinCluster::id( 0x22573 );

const double corSkinCluster::omega(DEFAULT_OMEGA);

const int corSkinCluster::_profileCategory(MProfiler::addCategory(NAME));

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

bool corSkinCluster::get_cor_valid(){
	return cor_valid;
}

void corSkinCluster::set_cor_valid(bool state){
	cor_valid = state;
}

void* corSkinCluster::creator()
{
	void *node = new corSkinCluster();
	((corSkinCluster *)node)->set_cor_valid(false);
	return node;
}

MStatus corSkinCluster::initialize()
{
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
				temp = weight_p[j]*weight_p[k]*weight_v[j]*weight_v[k];
				temp *= exp(-(pow(((weight_p[j]*weight_v[k]-weight_p[k]*weight_v[j])/pow(omega,2.0)),2.0)));
				result += temp;
			}
		}  // end k loop
	}  // end j loop

	return MStatus::kSuccess;
}

MStatus corSkinCluster::precomp(MDataBlock block)
{
	MStatus stat;

	// make sure the cor array is cleared out
	cor_ar.clear();

	MItMeshPolygon T(inputGeom, &stat);
	if (stat == MStatus::kFailure){
		std::cerr << "corSkinCluster::precomp, unable to get mesh iterator" << std::endl;
		return stat;
	}
	MItGeometry v_i(inputGeom, &stat);
	if (stat == MStatus::kFailure){
		std::cerr << "corSkinCluster::precomp, unable to get point iterator" << std::endl;
		return stat;
	}

	MArrayDataHandle w_i = block.inputArrayValue(weightList);
	if ( w_i.elementCount() == 0 ) {
		// no weights - nothing to do
		return MStatus::kFailure;
	}

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
		stat = T.numTriangles(num_tris);
		if (stat == MStatus::kSuccess){
			// for each triangle
			for (idx = 0; idx < num_tris; idx++){
				// get the verts
				stat = T.getTriangle(idx, tri_verts, tri_idx, MSpace::kObject);
				if (stat != MStatus::kSuccess){
					std::cerr << "corSkinCluster::precomp, unable to get triange" << std::endl;
					return MStatus::kFailure;
				}
				alpha = tri_verts[0];
				beta = tri_verts[1];
				gamma = tri_verts[2];
				// calc and store area of triange
				beta_alpha = MVector(beta-alpha);
				gamma_alpha = MVector(gamma-alpha);
				stat = tri_area.append(((beta_alpha ^ gamma_alpha).length())*0.5);
				if (stat != MStatus::kSuccess){
					std::cerr << "corskinCluster::precomp, unable to append area" << std::endl;
					return MStatus::kFailure;
				}
				
				// calc and store average vertex position
				stat = tri_avg_pos.append((alpha+beta+gamma)/3);
				if (stat != MStatus::kSuccess){
					std::cerr << "corSkinCluster::precomp, unable to apped average position" << std::endl;
				}
				
				// calc and store avg weights

				// get alpha weights
				stat = w_i.jumpToElement(tri_idx[0]);
				if (stat = MStatus::kSuccess){
					std::cerr << "corSkinCluster::precomp, unable to get weights for alpha: " << tri_idx[0] << std::endl;
					return MStatus::kFailure;
				}
				MArrayDataHandle alpha_weightsHandle = w_i.inputValue().child(weights);

				// get beta weights
				stat = w_i.jumpToElement(tri_idx[1]);
				if (stat != MStatus::kSuccess){
					std::cerr << "corSkinCluster::precomp, unable to get weights for alpha: " << tri_idx[1] << std::endl;
					return MStatus::kFailure;
				}
				MArrayDataHandle beta_weightsHandle = w_i.inputValue().child(weights);

				// get gamma weights
				stat = w_i.jumpToElement(tri_idx[2]);
				if (stat != MStatus::kSuccess){
					std::cerr << "corSkinCluster::precomp, unable to get weights for alpha: " << tri_idx[2] << std::endl;
					return MStatus::kFailure;
				}
				MArrayDataHandle gamma_weightsHandle = w_i.inputValue().child(weights);

				double a, b, c;
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
					tri_avg_weight[i] = (a + b + c)/3.0;
				} // end for
				tri_avg_weights.push_back(tri_avg_weight);
			} // end for
		}else{ // if num triangles fail
			std::cerr << "corSkinCluster::precomp, face has no triangles?" << std::endl;
			return MStatus::kFailure;
		} //end else
	} // end while

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
		MArrayDataHandle vertex_weights_handle = w_i.inputValue().child(weights);
		for (int i = 0; i < num_transforms; i++){
			stat = vertex_weights_handle.jumpToElement(i);
			if (stat != MStatus::kSuccess){
				vertex_weights[i] = 0.0;
			}else{
				vertex_weights[i] = vertex_weights_handle.inputValue().asDouble();
			}
		}

		// for each triangle
		for (int i = 0; i < num_tris; i++){
			s = 0.0;
			upper = MPoint(0.0,0.0,0.0);
			lower = 0.0;
			stat = similarity(vertex_weights, tri_avg_weights[i], num_transforms, s);
			upper += tri_avg_pos[i]*s*tri_area[i];
			lower += s*tri_area[i];
		}
		cor = upper/lower;
		cor_ar.append(cor);
		
		// iterate the loop
		v_i.next();
		w_i.next();
	} // end while

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
	if (!get_cor_valid()){
		precomp(block);
	}

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
	MQuaternion q;
	MQuaternion &rq = q;
	for (int j=0; j<numTransforms; ++j ) {
		const MMatrix &temp = transforms[j];
		rq = temp;
		q_j[j] = rq;
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

		MArrayDataHandle weightHandle = weightListHandle.inputValue().child(weights);

		MQuaternion q(0.0,0.0,0.0,0.0);
		MQuaternion result;
		double w_j;

		for (int j = 0; j < numTransforms;  j++){
			MStatus stat;
			stat = weightHandle.jumpToElement(j);
			if (stat == MStatus::kSuccess){
				w_j = weightHandle.inputValue().asDouble();
			}else{
				w_j = 0.0;
			}
			stat = qlerp(q,q_j[j], q);
		}

		q.normalizeIt();
		MMatrix R = q.asMatrix();

		//LBS Matrix
		MMatrix R_t_prime = MMatrix();
		for (int j = 0; j < numTransforms; j++){
			MTransformationMatrix tm(transforms[j]);
			MVector translate = tm.getTranslation(MSpace::kObject);
			MTransformationMatrix R_t = tm.asRotateMatrix();
			R_t.setTranslation(translate, MSpace::kObject);
			MMatrix temp = R_t.asMatrix();
			MStatus stat;
			stat = weightHandle.jumpToElement(j);
			if (stat == MStatus::kSuccess){
				w_j = weightHandle.inputValue().asDouble();
			}else{
				w_j = 0.0;
			}
			R_t_prime *= (temp * w_j);
		}

		MTransformationMatrix R_t_prime_tm(R_t_prime);

		MVector t = (cor_ar[i] * R_t_prime_tm.asRotateMatrix()) - (cor_ar[i] * R);
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
