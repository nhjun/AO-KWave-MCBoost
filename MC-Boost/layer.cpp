

#include "layer.h"
#include <iostream>
using std::cout;



Layer::Layer(Layer_Properties &props)
{
    this->mu_a = props.mu_a;
    this->mu_s = props.mu_s;
    this->mu_t = mu_a + mu_s;
    
    albedo = mu_s/(mu_s + mu_a);
    g = props.anisotropy;
    this->refractive_index = props.refractive_index;
    
    this->depth_start = props.start_depth;
    this->depth_end = props.end_depth;
}


Layer::Layer(double mu_a, double mu_s, double refractive_index, double anisotropy,
			 double depth_start, double depth_end)
{
	this->mu_a = mu_a;
	this->mu_s = mu_s;
	this->mu_t = mu_a + mu_s;
	albedo = mu_s/(mu_s + mu_a);
	g = anisotropy;
	this->refractive_index = refractive_index;
	
	this->depth_start = depth_start;
	this->depth_end = depth_end;
    
}

Layer::~Layer(void)
{
    // Free any memory allocated on the heap by this object.
    for (std::vector<Absorber *>::iterator i = p_absorbers.begin(); i < p_absorbers.end(); i++)
        delete *i;
}

void Layer::setAbsorpCoeff(double mu_a)
{
	this->mu_a = mu_a;
	
	// If we ever update the absorption coefficient we need to update the
	// transmission coefficient and albedo similarly.
	this->mu_t = mu_a + mu_s;
	updateAlbedo();
}

void Layer::setScatterCoeff(double mu_s)
{
	this->mu_s = mu_s;

	// If we ever update the scattering coefficient we need to update the
	// transmission coefficient albedo similarly.
	this->mu_t = mu_a + mu_s;
	updateAlbedo();
}


void Layer::updateAlbedo()
{
	albedo = mu_s/(mu_s + mu_a);
}


void Layer::addAbsorber(Absorber * absorber)
{
    // FIXME: Ensure the absorber fits within the bounds of the layer.
    //       
    p_absorbers.push_back(absorber);
}

	

// Returns the refractive index based on coordinates of photon.  That is,
// there could be an occlusion in the medium and we ensure the correct value
// is returned.
double Layer::getRefractiveIndex(const boost::shared_ptr<Vector3d> photonVector)
{
    cout << "Layer::getRefractiveIndex(Vector3d) stub\n";
}



// Returns the absorption coefficient after checking to see if the
// photon might be within an absorber.
double Layer::getAbsorpCoeff(const boost::shared_ptr<Vector3d> photonVector)
{
    // Iterate over all the absorbers in this layer and see if the coordinates
    // of the photon reside within the bounds of the absorber.  If so, we return
    // the absorption coefficient of the absorber, otherwise we return the 
    // absorption coefficient of the ambient layer.
    
    for (std::vector<Absorber *>::iterator it = p_absorbers.begin(); it != p_absorbers.end(); it++)
    {
        if ((*it)->inAbsorber(photonVector))
        {
            return (*it)->getAbsorberAbsorptionCoeff();
        }
    }
    
    // If we make it out of the loop (i.e. the photon is not in an absorber) we 
    // return the layer's absorption coefficient.
    return mu_a;
    
}


// Returns the absorption coefficient after checking to see if the
// photon might be within an absorber.
double Layer::getScatterCoeff(const boost::shared_ptr<Vector3d> photonVector)
{
    // Iterate over all the absorbers in this layer and see if the coordinates
    // of the photon reside within the bounds of the absorber.  If so, we return
    // the absorption coefficient of the absorber, otherwise we return the
    // absorption coefficient of the ambient layer.

    for (std::vector<Absorber *>::iterator it = p_absorbers.begin(); it != p_absorbers.end(); it++)
    {
        if ((*it)->inAbsorber(photonVector))
        {
            return (*it)->getAbsorberScatteringCoeff();
        }
    }

    // If we make it out of the loop (i.e. the photon is not in an absorber) we
    // return the layer's absorption coefficient.
    return mu_s;

}


double Layer::getTotalAttenuationCoeff(const boost::shared_ptr<Vector3d> photonVector)
{
    // Iterate over all the absorbers in this layer and see if the coordinates
    // of the photon reside within the bounds of the absorber.  If so, we return
    // the total attenuation coefficient of the absorber, otherwise we return the
    // total attenuation coefficient of the ambient layer.
    
    for (std::vector<Absorber *>::iterator it = p_absorbers.begin(); it != p_absorbers.end(); it++)
    {
        if ((*it)->inAbsorber(photonVector))
        {
            return ((*it)->getAbsorberScatteringCoeff() + (*it)->getAbsorberScatteringCoeff());
        }
    }
    
    // If we make it out of the loop (i.e. the photon is not in an absorber) we
    // return the layer's total attenuation coefficient.
    return (mu_a + mu_s);
}


void Layer::updateAbsorbedWeightByAbsorber(const boost::shared_ptr<Vector3d> photonVector, const double absorbed)
{
    // Iterate over all the absorbers in this layer and see if the coordinates
    // of the photon reside within the bounds of the absorber.  If so, we return
    // the absorption coefficient of the absorber, otherwise we return the 
    // absorption coefficient of the ambient layer.
    
    for (std::vector<Absorber *>::iterator it = p_absorbers.begin(); it != p_absorbers.end(); it++)
    {
        if ((*it)->inAbsorber(photonVector))
        {
            (*it)->updateAbsorbedWeight(absorbed);
        }
    }
}

Absorber * Layer::getAbsorber(const boost::shared_ptr<Vector3d> photonVector)
{
    // Iterate over all the absorbers in this layer and see if the coordinates
    // of the photon reside within the bounds of the absorber.  If so, we return
    // the absorption coefficient of the absorber, otherwise we return the 
    // absorption coefficient of the ambient layer.
    for (std::vector<Absorber *>::iterator it = p_absorbers.begin(); it != p_absorbers.end(); it++)
    {
        if ((*it)->inAbsorber(photonVector))
        {
            return *it;
        }
    }
    
    return NULL;
}

// Iterate over all absorbers and write their data out to file.
void Layer::writeAbsorberData(void)
{
    // Write out the data for every absorber in the medium.
    for (std::vector<Absorber *>::iterator it = p_absorbers.begin(); it != p_absorbers.end(); it++)
    {
        (*it)->writeData();
        
    }
}

