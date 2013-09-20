// Defines attributes of a layer.
#ifndef LAYER_H
#define LAYER_H


#include "absorber.h"
#include "coordinates.h"
#include <vector>



typedef struct {
    
    double mu_a;
	double mu_s;
	double refractive_index;
	double anisotropy;
	double start_depth;
	double end_depth;
    
} Layer_Properties;



class Layer
{	

public:
    Layer(Layer_Properties &properties);
	Layer(double mu_a, double mu_s, double ref_index, double anisotropy,
		  double depth_start, double depth_end);
	~Layer(void);


    // Returns the absorption coefficient of the layer.
	double	getAbsorpCoeff(void) const	{return mu_a;}
    // Returns the absorption coeffiecient of the layer based on the photon's coordinates
    // Checks are made to see if the photon has made it's way into an absorber as well.
    double  getAbsorpCoeff(const boost::shared_ptr<Vector3d> photonVector);
    
    // Returns the scattering coefficient of the layer.
	double	getScatterCoeff(void) const	{return mu_s;}
    double  getScatterCoeff(const boost::shared_ptr<Vector3d> photonVector);
    
    // Returns total interaction coefficient (mu_a + mu_s).
	double	getTotalAttenuationCoeff(void) const	{return mu_t;}
    double  getTotalAttenuationCoeff(const boost::shared_ptr<Vector3d> photonVector);
    
    // Return the albedo
	double	getAlbedo(void) const			{return albedo;}
    
	// Return the anisotropy of the layer.
	double	getAnisotropy(void) 		{return g;}
    double  getAnisotropy(const boost::shared_ptr<Vector3d> photonVector);
    

    // Return the impedance of the layer.
    double 	getImpedance(void) {return impedance;}

	double 	getDepthStart(void) const 		{return depth_start;}
	double  getDepthEnd(void)	const		{return depth_end;}
    
	// Return the refractive index of the layer.
	double	getRefractiveIndex(void) const	{return refractive_index;}
    double  getRefractiveIndex(const boost::shared_ptr<Vector3d> photonVector);

	void	setAbsorpCoeff(const double mu_a);
	void	setScatterCoeff(const double mu_s);

    /// Any time the optical properties of the medium are set or changed, this needs to be called.
    /// to propagate those changes to all params (ie. albedo, total attenuation, etc.
    void	Update_optical_params();
    
    void    addAbsorber(Absorber * a);
    
    void    updateAbsorbedWeightByAbsorber(const boost::shared_ptr<Vector3d> currLocation, const double absorbed);
    
    // Iterate over all absorbers and write their data out to file.
    void    writeAbsorberData(void);
    
    // Return the absorber at this location 'currLocation' in the medium.
    Absorber * getAbsorber(const boost::shared_ptr<Vector3d> currLocation);

    /// Return the number of absorbers in a layer.
    size_t  Get_num_absorbers(void) {return p_absorbers.size();};
    

	
private:
    
	// Anisotropy factor.
	double g;
	
	// Absorption coefficient
	double mu_a;
	
	// Scattering coefficient
	double mu_s;
	
	// Transmission coefficient
	double mu_t;
	
	// The refractive index of the layer
	double refractive_index;
	
	// The width of the layer.
	//double radial_size;
	
	// z-coordinate value at which the layer starts.
	double depth_start;
	
	// z-coordinate value at which the layer ends.
	double depth_end;
	
	// Albedo of the layer.
	double albedo;	

	// The impedance of the layer.
	double impedance;
    
    // A vector that holds all the abosrbers in this layer.
    std::vector<Absorber *> p_absorbers;
	
};

#endif // end LAYER_H


