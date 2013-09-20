
#include "photon.h" // Photon class is a friend of the Medium class.
#include "debug.h"
#include "vector3D.h"
#include "layer.h"
#include "medium.h"
#include "detector.h"
#include "pressureMap.h"
#include "refractiveMap.h"
#include "displacementMap.h"
#include <cmath>
#include <cassert>
#include <cstdlib>

#undef DEBUG

Medium::Medium()
{
	z_bound = x_bound = y_bound = 0.0f; // Default bounds of the medium [meters].
	this->initCommon();
}

Medium::Medium(const double x, const double y, const double z)
{
	this->z_bound = z;
	this->y_bound = y;
	this->x_bound = x;
	this->initCommon();
}

Medium::~Medium()
{
	// Free the PressureMap() object.
	if (kwave.pmap) {
		delete kwave.pmap;
		kwave.pmap = NULL;
	}

	// Free the RefractiveMap() object.
	if (kwave.nmap) {
		delete kwave.nmap;
		kwave.nmap = NULL;
	}

	// Free the DisplacementMap() object.
	if (kwave.dmap) {
		delete kwave.dmap;
		kwave.dmap = NULL;
	}

	// If there were any absorbers in the medium, write out their data.
	for (vector<Layer *>::iterator it = p_layers.begin(); it != p_layers.end(); it++)
	{
		(*it)->writeAbsorberData();
		delete *it;
	}
}


void Medium::initCommon(void)
{

    kwave.pmap = NULL;     	// Pointer to a PressureMap() object.
	kwave.nmap = NULL;		// Pointer to a RefractiveMap() object.
	kwave.dmap = NULL;		// Pointer to a DisplacementMap() object.

    X_PML_OFFSET = 0;
    Y_PML_OFFSET = 0;
    Z_PML_OFFSET = 0;

    //coords_file.open("photon-paths.txt");
	//photon_data_file.open("photon-exit-data.txt");

}


// Add the layer to the medium by pushing it onto the vector container.
void Medium::addLayer(Layer *layer)
{
    /// Display the properties of the added layer in the medium.
    cout << "-----------------------------------------------------\n"
         << "Adding a layer to the medium /\n"
         << "-----------------------------\n"
         << " absorption = " << layer->getAbsorpCoeff()/100 << " [cm^-1]\n"
         << " scattering = " << layer->getScatterCoeff()/100 << " [cm^-1]\n"
         << " anisotropy = " << layer->getAnisotropy() << '\n'
         << " refractive index = " << layer->getRefractiveIndex() << '\n'
         << " number of absorbers = " << layer->Get_num_absorbers() << '\n';

	p_layers.push_back(layer);
}
void Medium::addLayer(Layer_Properties props)
{
    /// Display the properties of the added layer in the medium.
    cout << "-----------------------------------------------------\n"
         << "Adding a layer to the medium /\n"
         << "-----------------------------\n"
         << " absorption = " << props.mu_a/100 << " [cm^-1]\n"
         << " scattering = " << props.mu_s/100 << " [cm^-1]\n"
         << " anisotropy = " << props.anisotropy << '\n'
         << " refractive index = " << props.refractive_index << '\n';

    p_layers.push_back(new Layer(props));
}




// Add an 3-dimensional pressure map object to the medium.
void Medium::addPressureMap(PressureMap *p_map)
{
	assert(p_map != NULL);

    /// Since we are adding a new pressure map, there is a chance one is already assigned.
    /// If one already exists, we free the memory and assign the new one.
    if (kwave.pmap != NULL)
    {
        delete kwave.pmap;
        kwave.pmap = NULL;
    }

    /// Assign new pressure map.
	kwave.pmap = p_map;

}

// Add a 3D refractive map object to the medium.
void Medium::addRefractiveMap(RefractiveMap *n_map)
{
	assert(n_map != NULL);

    /// Since we are adding a new refractive map, there is a chance one is already assigned.
    /// If one already exists, we free the memory and assign the new one.
    if (kwave.nmap != NULL)
    {
        delete kwave.nmap;
        kwave.nmap = NULL;
    }

    /// Assign new refractive map.
	kwave.nmap = n_map;

}


// Add an 3-dimensional displacement map object to the medium.
void Medium::addDisplacementMap(DisplacementMap *d_map)
{
	assert(d_map != NULL);

    /// Since we are adding a new displacement map, there is a chance one is already assigned.
    /// If one already exists, we free the memory and assign the new one.
    if (kwave.dmap != NULL)
    {
        delete kwave.dmap;
        kwave.dmap = NULL;
    }

    /// Assign the new displacement map.
	kwave.dmap = d_map;
}





void Medium::addDetector(Detector *detector)
{

	p_detectors.push_back(detector);
}



// See if photon has crossed the detector plane.
int Medium::photonHitDetectorPlane(const boost::shared_ptr<Vector3d> p0)
{
	bool hitDetectorNumTimes = 0;
	// Free the memory for layers that were added to the medium.
	for (vector<Detector *>::iterator it = p_detectors.begin(); it != p_detectors.end(); it++)
	{
		if ((*it)->photonHitDetector(p0))
			hitDetectorNumTimes++;
	}

	return hitDetectorNumTimes;
}

Layer * Medium::getLayerAboveCurrent(Layer *currentLayer)
{
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(currentLayer != NULL);

	// If we have only one layer, no need to iterate through the vector.
	// And we should return NULL since there is no layer above us.
	if (p_layers.size() == 1)
		return NULL;



	// Otherwise we walk the vector and return 'trailer' since it is the
	// one before the current layer (i.e. 'it').
	vector<Layer *>::iterator it;
	vector<Layer *>::iterator trailer;
	it = p_layers.begin(); // Get the first layer from the array.

	// If we are at the top of the medium there is no layer above, so return NULL;
	if (currentLayer == (*it))
		return NULL;


	///boost::mutex::scoped_lock lock(m_layer_above_mutex);
	while(it != p_layers.end()) {
		trailer = it;  // Assign the trailer to the current layer.
		it++;         // Advance the iterator to the next layer.

		// Find the layer we are in within the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer
		// because trailer will be pointing to the previous layer in the medium.
		//if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z)
		if ((*it) == currentLayer)
			break;
	}

	// Sanity check.  If the trailer has made it to the end, which means
	// the iterator made it past the end, then there
	// was no previous layer found, and something went wrong.
	if (trailer == p_layers.end())
		return NULL;


	// If we make it here, we have found the previous layer.
	return *trailer;
}


Layer * Medium::getLayerBelowCurrent(double z)
{
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(z >= 0 && z <= z_bound);

	// If we have only one layer, no need to iterate through the vector.
	// And we should return NULL since there is no layer below us.
	if (p_layers.size() == 1)
		return NULL;

	// The case where there is no layer below is since we are at the bottom of the
	// medium.
	if (z == z_bound)
		return NULL;



	vector<Layer *>::iterator it;
    ///boost::mutex::scoped_lock lock(m_layer_below_mutex);
	for (it = p_layers.begin(); it != p_layers.end(); it++) {
		// Find the layer we are in within the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer.
		if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z) {
			return *(++it);
		}
	}

	// If the above loop never returned a layer it means we made it through the list
	// so there is no layer below us, therefore we return null.
	return NULL;


}


// Return the layer in the medium at the passed in depth 'z'.
// We iterate through the vector which contains pointers to the layers.
// When the correct layer is found from the depth we return the layer object.
Layer * Medium::getLayerFromDepth(double z)
{
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(z >= 0 && z <= z_bound);


	vector<Layer *>::iterator it;
	for (it = p_layers.begin(); it != p_layers.end(); it++) {
		// Find the layer we are in within the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer.
		if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z)
			break;
	}

	// Return layer based on the depth passed in.
	return *it;
}


double Medium::getLayerAbsorptionCoeff(double z)
{
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(z >= 0 && z <= z_bound);

	double absorp_coeff = -1;
	vector<Layer *>::iterator it;

	for (it = p_layers.begin(); it != p_layers.end(); it++) {
		// Find the layer we are it in the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer.
		if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z) {
			absorp_coeff = (*it)->getAbsorpCoeff();
			break;
		}
	}

	// If not found, report error.
	assert(absorp_coeff != 0);

	// If not found, fail.
	// If not found, report error.
	assert(absorp_coeff != -1);

	// Return the absorption coefficient value.
	return absorp_coeff;
}


double Medium::getLayerScatterCoeff(double z)
{
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(z >= 0 && z <= z_bound);

	double scatter_coeff = -1;
	vector<Layer *>::iterator it;


	for (it = p_layers.begin(); it != p_layers.end(); it++) {
		// Find the layer we are it in the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer.
		if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z) {
			scatter_coeff = (*it)->getScatterCoeff();
			break;
		}
	}

	// If not found, report error.
	assert(scatter_coeff != 0);

	// If not found, fail.
	// If not found, report error.
	assert(scatter_coeff != -1);

	// Return the scattering coefficient for the layer that resides at depth 'z'.
	return scatter_coeff;
}


double Medium::getAnisotropyFromDepth(double z)
{
	// Ensure that the photon's z-axis coordinate is sane.  That is,
	// it has not left the medium.
	assert(z >= 0 && z <= z_bound);

	double anisotropy = -1;
	vector<Layer *>::iterator it;


	for (it = p_layers.begin(); it != p_layers.end(); it++) {
		// Find the layer we are it in the medium based on the depth (i.e. z)
		// that was passed in.  Break from the loop when we find the correct layer.
		if ((*it)->getDepthStart() <= z && (*it)->getDepthEnd() >= z) {
			anisotropy = (*it)->getAnisotropy();
			break;
		}
	}

	// If not found, report error.
	assert(anisotropy != 0);

	// If not found, fail.
	// If not found, report error.
	assert(anisotropy != -1);

	// Return the anisotropy value for the layer that resides at depth 'z'.
	return anisotropy;
}





/// Create the displacement map based on what is returned from the computation
/// that happens in KSpaceSolver.
void Medium::Create_displacement_map(TRealMatrix * disp_x,
                                     TRealMatrix * disp_y,
                                     TRealMatrix * disp_z)
{



    if (kwave.dmap == NULL)
    {

        kwave.dmap = new DisplacementMap();
        kwave.dmap->Assign_displacement_map(disp_x,
                                            disp_y,
                                            disp_z);
    }
    else
    {
        kwave.dmap->Assign_displacement_map(disp_x,
                                            disp_y,
                                            disp_z);
    }


}



/// Create the refractive index map based on what is returned from the computation
/// that happens in KSpaceSolver or from what was loaded in from an HDF5 file from a previous run.
void Medium::Create_refractive_map(TRealMatrix * refractive_total)
{
    /// XXX:
    /// - This is defined here, but would change depending on temperature changes, if that
    ///   is ever incorporated into this simulation, which in turn would affect density and SOS.
    background_refractive_index = 1.33;

    if(kwave.nmap == NULL)
    {
        /// Takes pressure as an argument in the constructor and forms the refractive map data.
        kwave.nmap = new RefractiveMap();
        kwave.nmap->Assign_refractive_map(refractive_total);
    }
    else
    {
        /// Updates the refractive map data.
        kwave.nmap->Assign_refractive_map(refractive_total);
    }

}


void Medium::Create_refractive_map(TRealMatrix * refractive_x,
                                   TRealMatrix * refractive_y,
                                   TRealMatrix * refractive_z)
{

    /// XXX:
    /// - This is defined here, but would change depending on temperature changes, if that
    ///   is ever incorporated into this simulation, which in turn would affect density and SOS.
    background_refractive_index = 1.33;


    if(kwave.nmap == NULL)
    {
        /// Takes pressure as an argument in the constructor and forms the refractive map data.
        kwave.nmap = new RefractiveMap();
        kwave.nmap->Assign_refractive_map(refractive_x,
                                          refractive_y,
                                          refractive_z);
    }
    else
    {
        /// Updates the refractive map data.
        kwave.nmap->Assign_refractive_map(refractive_x,
                                          refractive_y,
                                          refractive_z);
    }

}


