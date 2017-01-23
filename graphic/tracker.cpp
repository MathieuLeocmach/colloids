/**
    Copyright 2008,2009 Mathieu Leocmach

    This file is part of Colloids.

    Colloids is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Colloids is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Colloids.  If not, see <http://www.gnu.org/licenses/>.

 * \file tracker.cpp
 * \brief Implement classes to track particles in 3D image series
 * \author Mathieu Leocmach
 * \date 7 april 2009
 *
 * Object oriented implementation of the code elaborated in 2008.
 * Old auto-configuration routine are commented at the end but not implemented this time
 * Based on Crocker & Grier algorithm
 *
 */

#include <list>
#include <boost/progress.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include "tracker.hpp"
#include "particles.hpp"

#ifndef TRACKER_N_THREADS
#define TRACKER_N_THREADS 1
#endif

using namespace std;
using namespace Colloids;
//using namespace cimg_library;
typedef boost::multi_array_types::index_range range;

/** @brief Tracker Constructor  */
 Tracker::Tracker(const boost::array<size_t,3> &dims, const unsigned fs)
{
	setFlags(fs);
	setDimensions(dims);
}

/** @brief Destructor  */
 Tracker::~Tracker()
{
	fftwf_destroy_plan(forward_plan);
	fftwf_destroy_plan(backward_plan);
	fftwf_free(data);
}

/** @brief Set the dimensions of the image, allocate memory, etc.
	\param dims Sizes of each non trival dimension (min=1, max=3), supposed in row major order.
	n[0] varies slower than n[1], itself variing slower that n[3]
	For a 3D image scanned with x faster than y faster than z, the dimensions must be given in reverse order
	n[0]=dimz, n[1]=dimy, n[2]=dimx
*/
void Tracker::setDimensions(const boost::array<size_t,3> &dims)
{
	//allocate main memory block for FFT.
	// Last dimension has to be padded with extra values to allow real2complex and c2r fft
	boost::array<size_t,3> paddedDims = dims;
	paddedDims.back()= 2*(paddedDims.back()/2+1);
	size_t memsize = 1;
	memsize = accumulate(paddedDims.begin(),paddedDims.end(),1,multiplies<size_t>());
	cout<<"Allocating a block of "<<sizeof(float) * memsize<<" bytes ... ";
	data = (float*)fftwf_malloc(sizeof(float)* memsize);
	assert(data);

	//allocate memory.
	centersMap.resize(dims);
	paddedDims.back() = dims.back()/2 + 1;
	FFTmask.resize(paddedDims);

	//planning fft.
	int n[3];
	copy(dims.begin(),dims.end(),&(n[0]));
	forward_plan = fftwf_plan_dft_r2c(dims.size(), &(n[0]), data, (fftwf_complex *)data, flags);
	//n[2] = 2*(n[2]/2 +1);
	backward_plan = fftwf_plan_dft_c2r(dims.size(), &(n[0]), (fftwf_complex *)data, data, flags);
}

/** @brief make the band pass filter mask in Fourier space
    \param radiiMin Real space mimimum radiii in each dimension, supposed in row major order
    \param radiiMax Real space maximum radiii in each dimension, supposed in row major order
*/
void Tracker::makeBandPassMask(const boost::array<double,3> &radiiMin, const boost::array<double,3> &radiiMax)
{
    if(!quiet)
    {
    	cout << "making the band-pass filter (";
		copy(FFTmask.shape(),FFTmask.shape()+FFTmask.num_dimensions(),ostream_iterator<size_t>(cout,","));
		cout<<") ... ";
    }

    if(!equal(radiiMin.begin(),radiiMin.end(),radiiMax.begin(),less<double>()))
		throw invalid_argument("For each dimension we must have radiusMin < radiusMax ");

    boost::array<double,3> freqMax, freqMin;
    transform(centersMap.shape(),centersMap.shape()+centersMap.num_dimensions(),radiiMin.begin(), freqMax.begin(), divides<double>());
    transform(centersMap.shape(),centersMap.shape()+centersMap.num_dimensions(),radiiMax.begin(), freqMin.begin(), divides<double>());
    transform(freqMax.begin(), freqMax.end(), freqMax.begin(), bind2nd(divides<double>(), 2.0));
    transform(freqMin.begin(), freqMin.end(), freqMin.begin(), bind2nd(divides<double>(), 2.0));
    if(!quiet)
    {
        cout << "Size min (";
		copy(radiiMin.begin(),radiiMin.end(),ostream_iterator<double>(cout,","));
		cout<<") ... ";
		cout << "Size max (";
		copy(radiiMax.begin(),radiiMax.end(),ostream_iterator<double>(cout,","));
		cout<<") ... ";
    	cout << "Freq max (";
		copy(freqMax.begin(),freqMax.end(),ostream_iterator<double>(cout,","));
		cout<<") ... ";
		cout << "Freq min (";
		copy(freqMin.begin(),freqMin.end(),ostream_iterator<double>(cout,","));
		cout<<") ... ";
    }

    //The spectrum output of FFT is centered on (0,0,0), not in the center of the image.
    //The last dimension is only (n/2+1), with n the last dimension of the real-space image
    //So the band pass mask is half on a 3D-ellipsoidal shell centered on (0,0,0).
    //Moreover, the symetry planes allow to draw only 1/2^d th of the mask and then mirror it.

    //get a view of 1/2^d th of the mask

    boost::array<size_t,3> halfdims;
    for(int d=0;d<3;++d)
		halfdims[d] = centersMap.shape()[d] / 2 + 1;
    array_type_b::array_view<3>::type portion = FFTmask[
		boost::indices[range(0,halfdims[0])][range(0,halfdims[1])][range(0,halfdims[2])]
		];
	if(!quiet)
    {
    	cout<<"portion (";
		copy(portion.shape(),portion.shape()+portion.num_dimensions(),ostream_iterator<size_t>(cout,","));
		cout<<") ... ";
    }
    //boost::multi_array_ref<bool, 3>portion(data, halfdims);
    double imin, imax, jmin, jmax;
    size_t di=0,dj=0;
    //array_type_b::array_view<3>::type::array_view<1>::type line;

	//draw a 1/2^d th of the mask
	//not the best efficiency if z is flat, but then the number of pixels is small so no problem
	for(array_type_b::iterator i = portion.begin();i!=portion.end();++i)
	{
		dj=0;
		imin = pow(di / freqMin[0], 2.0);
		imax = pow(di / freqMax[0], 2.0);
		for(array_type_b::value_type::iterator j = i->begin();j!=i->end();++j)
		{
		    //cout<<(j-i->begin())<<endl;
		    jmin = pow(dj / freqMin[1], 2.0);
		    jmax = pow(dj / freqMax[1], 2.0);
			size_t kmin = min(
				j->size(),
				(imin+jmin>1.0)?0:
				(size_t)(1+freqMin[2]*sqrt(1.0 - imin - jmin))
				);
			size_t kmax = min(
				j->size(),
				(imax+jmax>1.0)?kmin:
				(size_t)(1+freqMax[2]*sqrt(1.0 - imax - jmax))
				);
			//cout<<imin<<","<<imax<<","<<jmin<<","<<jmax<<","<<kmin<<","<<kmax;//<<endl;
			//cout<<" 0->"<<kmin;
			fill_n(j->begin(), kmin, false);
			//cout<<"->"<<kmax;
			fill(j->begin() + kmin, j->begin() + kmax, true);
			//cout<<"->"<<j->size()<<endl;
			fill(j->begin() + kmax, j->end(), false);
			dj++;
		}
		di++;
	}

	this->mirrorMask();

	if(view) displayMask();
}

/** @brief make a low pass filter mask in Fourier space (to remove high frequency noise
    \param radiiMin Real space mimimum radiii in each dimension, supposed in row major order
*/
void Tracker::makeLowPassMask(const boost::array<double,3> &radiiMin)
{
    if(!quiet)
    {
    	cout << "making the high-pass filter (";
		copy(FFTmask.shape(),FFTmask.shape()+FFTmask.num_dimensions(),ostream_iterator<size_t>(cout,","));
		cout<<") ... ";
    }

    boost::array<double,3> freqMax;
    transform(centersMap.shape(),centersMap.shape()+centersMap.num_dimensions(),radiiMin.begin(), freqMax.begin(), divides<double>());
    transform(freqMax.begin(), freqMax.end(), freqMax.begin(), bind2nd(divides<double>(), 2.0));
    if(!quiet)
    {
        cout << "Size min (";
		copy(radiiMin.begin(),radiiMin.end(),ostream_iterator<double>(cout,","));
		cout<<") ... ";
    	cout << "Freq max (";
		copy(freqMax.begin(),freqMax.end(),ostream_iterator<double>(cout,","));
		cout<<") ... ";
    }

    //The spectrum output of FFT is centered on (0,0,0), not in the center of the image.
    //The last dimension is only (n/2+1), with n the last dimension of the real-space image
    //So the band pass mask is half on a 3D-ellipsoidal shell centered on (0,0,0).
    //Moreover, the symetry planes allow to draw only 1/2^d th of the mask and then mirror it.

    //get a view of 1/2^d th of the mask

    boost::array<size_t,3> halfdims;
    for(int d=0;d<3;++d)
		halfdims[d] = centersMap.shape()[d] / 2 + 1;
    array_type_b::array_view<3>::type portion = FFTmask[
		boost::indices[range(0,halfdims[0])][range(0,halfdims[1])][range(0,halfdims[2])]
		];
	if(!quiet)
    {
    	cout<<"portion (";
		copy(portion.shape(),portion.shape()+portion.num_dimensions(),ostream_iterator<size_t>(cout,","));
		cout<<") ... ";
    }
    //boost::multi_array_ref<bool, 3>portion(data, halfdims);
    double imax, jmax;
    size_t di=0,dj=0;
    //array_type_b::array_view<3>::type::array_view<1>::type line;

	//draw a 1/2^d th of the mask
	//not the best efficiency if z is flat, but then the number of pixels is small so no problem
	for(array_type_b::iterator i = portion.begin();i!=portion.end();++i)
	{
		dj=0;
		imax = pow(di / freqMax[0], 2.0);
		for(array_type_b::value_type::iterator j = i->begin();j!=i->end();++j)
		{
		    //cout<<(j-i->begin())<<endl;
		    jmax = pow(dj / freqMax[1], 2.0);
			size_t kmax = min(
				j->size(),
				(imax+jmax>1.0)?0:
				(size_t)(1+freqMax[2]*sqrt(1.0 - imax - jmax))
				);
			//cout<<" 0->"<<kmax;
			fill_n(j->begin(), kmax, true);
			//cout<<"->"<<j->size()<<endl;
			fill(j->begin() + kmax, j->end(), false);
			dj++;
		}
		di++;
	}

	this->mirrorMask();


	if(view) displayMask();
}

/** @brief mirror the first 1/2^d th of the FFT Mask in order to form the full mask  */
void Tracker::mirrorMask()
{
    boost::array<size_t,3> halfdims;
    for(int d=0;d<3;++d)
		halfdims[d] = centersMap.shape()[d] / 2 + 1;

    //mirror in first direction
	if(halfdims[0] >=3)
	{
        size_t parity = 1-(int)(FFTmask.shape()[0] % 2);
        array_type_b::array_view<3>::type original = FFTmask[
            boost::indices[range(parity, halfdims[0] -1 + parity)][range(0, halfdims[1])][range()]
            ];
        array_type_b::array_view<3>::type mirror = FFTmask[
            boost::indices[range(FFTmask.shape()[0]-1, FFTmask.shape()[0]-1-original.shape()[0], -1)][range(0, halfdims[1])][range()]
            ];
        if(!quiet)
        {
            cout<<"original(";
            copy(original.shape(),original.shape()+original.num_dimensions(),ostream_iterator<size_t>(cout,","));
            cout<<") -> mirror(";
            copy(mirror.shape(),mirror.shape()+mirror.num_dimensions(),ostream_iterator<size_t>(cout,","));
            cout<<") ... ";
        }
        mirror = original;
	}

	//mirror in second direction
	if(halfdims[1] >=3)
	{
        size_t parity =  1-(int)(FFTmask.shape()[1] % 2);
        array_type_b::array_view<3>::type original2 = FFTmask[
            boost::indices[range()][range(parity, halfdims[1] -1 + parity)][range()]
            ];
        array_type_b::array_view<3>::type mirror2 = FFTmask[
            boost::indices[range()][range(FFTmask.shape()[1]-1, FFTmask.shape()[1]-1-original2.shape()[1], -1)][range()]
            ];
        if(!quiet)
        {
            cout<<"original2(";
            copy(original2.shape(),original2.shape()+original2.num_dimensions(),ostream_iterator<size_t>(cout,","));
            cout<<") -> mirror2(";
            copy(mirror2.shape(),mirror2.shape()+mirror2.num_dimensions(),ostream_iterator<size_t>(cout,","));
            cout<<") ... ";
        }
        mirror2 = original2;
	}
}



/** @brief Get a view of the main data array, including padding  */
Tracker::array_type_r Tracker::get_padded_image(const boost::array<size_t,3>& ordering)
{
	boost::array<size_t,3> paddedExt;
	for(size_t d=0;d<3;++d)
		paddedExt[ordering[d]] = (d==2) ? 2*FFTmask.shape()[d] : centersMap.shape()[d];

	bool ascending[] = {true,true,true};
	return array_type_r(data, paddedExt, boost::general_storage_order<3>(ordering.begin(),ascending));
}

/** @brief get_image Get a view of the main data array  */
Tracker::view_type Tracker::get_image(const boost::array<size_t,3>& ordering)
{
	boost::array<size_t,3> paddedExt;
	boost::array<range,3> ranges;
	for(size_t d=0;d<3;++d)
	{
		ranges[ordering[d]] = range(0, centersMap.shape()[d]);
		paddedExt[ordering[d]] = (d==2) ? 2*FFTmask.shape()[d] : ranges[ordering[d]].finish();
	}
	/*cout<<"view (";
	for(size_t d=0;d<3;++d)
		cout<<ranges[d].finish()<<",";
	cout<<") of (";
	copy(paddedExt.begin(),paddedExt.end(),ostream_iterator<size_t>(cout,","));
	cout<<")"<<endl;*/

	bool ascending[] = {true,true,true};
	return array_type_r(data, paddedExt, boost::general_storage_order<3>(ordering.begin(),ascending))[
		boost::indices[ranges[0]][ranges[1]][ranges[2]]];
}

/** @brief unpad the last dimension of the image after inplace FFT  */
void Tracker::unpad()
{
	if(this->padded)
	{
		copyImage(data);
		this->padded = false;
	}
}

/** @brief normalize
  *
  * @todo: document this function
  */
void Tracker::normalize()
{
    boost::progress_timer ptimer;
    if(!quiet) cout << "normalize ... ";
    if(this->padded)
    {
        //padding pixels are not initialized. At best they are 0, but we can't count on it.
        //So a max_element of min_element would have no meaning.
        //We just divide by the number of pixels (with FFTW, ifft(fft(x)) = N*x)
        transform(data, data+2*FFTmask.num_elements(), data, bind2nd(divides<float>(), (float)centersMap.num_elements()*2));
        transform(data, data+2*FFTmask.num_elements(), data, bind2nd(plus<float>(),(float)128));
    }
    else
    {
        /*using namespace boost::accumulators;
        accumulator_set<long double, features<tag::min, tag::max, tag::mean> > acc;
        acc = for_each(data, data+centersMap.num_elements(), acc);
        const float minimum = boost::accumulators::min(acc), maximum = boost::accumulators::max(acc)*/
        std::pair<float*, float*> minmax = boost::minmax_element(data, data+centersMap.num_elements());
        /*if(!quiet) cout << "min="<< minimum
                        <<" max="<< maximum
                        //<<" mean="<< boost::accumulators::mean(acc)
                        <<" ... ";*/
        transform(data, data+centersMap.num_elements(), data, bind2nd(minus<float>(), *minmax.first));
        transform(data, data+centersMap.num_elements(), data, bind2nd(multiplies<float>(), 255.0/(*minmax.second - *minmax.first)));
    }
}

/** @brief set a reference brightness in order to adapt the threshold  */
/*void Tracker::setRefBrightness()
{
    if(this->padded)
        this->refBrightness = accumulate(
            data, data+2*FFTmask.num_elements(),
            0.0, plus<double>()
            ) * 2*FFTmask.num_elements() / (double)centersMap.num_elements();
    else
        this->refBrightness = accumulate(
            data, data+centersMap.num_elements(),
            0.0, plus<double>()
            );
    this->hasRefBrightness = true;
}*/



/** @brief Extract the coordinates of the local maxima of the filtered image
  *
  * Use :
	- After making the FFTmask (once) and loading the image data (every time step)
	- Before scaling the coordinates by Zratio
  */
Particles Tracker::trackXYZ(const float &threshold, bool autoThreshold)
{
    if(view) display("original stack");
	if(!quiet) cout << "band-passing ... ";
    FFTapplyMask();
    if(!quiet) cout << "done!" << endl;
    if(view) display("band passed");

    if(!quiet) cout << "pixel centers ... ";
    findPixelCenters(threshold);
    if(view) displayCenters();

    if(!quiet) cout << "sub-pixel resolution ... ";
    return getSubPixel(autoThreshold);
}

/** @brief get the central pixel brightess for each coordinate in centers.  */
vector<float> Tracker::getIntensities(const Particles &centers)
{
    vector<float> ret(centers.size());
    for(size_t p=0; p<centers.size(); ++p)
    {
        //convert real coordinates to integer strides
        boost::array<size_t, 3> tmp;
        for(size_t d=0; d<tmp.size(); ++d)
            tmp[d] = centersMap.strides()[d] * (size_t)centers[p][(fortran_order)?(2-d):d];

        //normlize the value by the size of the image to get rid of the FTTW unnormalization and be in the same unit as the threshold parameter
        ret[p] = (*(data+accumulate(tmp.begin(), tmp.end(), 0)))/centersMap.num_elements();
    }
    return ret;
}



/** @brief load FFTW wisdom from file
  * \return	false if the file isn't found
  */
bool Tracker::loadWisdom(const std::string & filename)
{
	FILE * wis = fopen(filename.c_str(),"r");
	if(wis == NULL) return false;
	fftwf_import_wisdom_from_file(wis);
	fclose(wis);
	return true;
}

/** @brief save FFTW wisdom (plans for fastest FFT)  */
void Tracker::saveWisdom(const std::string & filename)
{
	//no overwritting of the content of an existing file
	loadWisdom(filename);
	FILE * wis = fopen(filename.c_str(),"w+");
	if(wis == NULL)
		throw std::invalid_argument("No such file as "+filename);
	fftwf_export_wisdom_to_file(wis);
	fclose(wis);
}

/** @brief export the mask to file as ascii  */
void Tracker::maskToFile(std::string &filename) const
{
	ofstream os(filename.c_str());
	copy(
		FFTmask.origin(),FFTmask.origin()+FFTmask.num_elements(),
		 ostream_iterator<bool>(os,	"\t")
		);
	os.close();
}

/** @brief export the image to file as binary  */
void Tracker::imageToFile(std::string &filename) const
{
	ofstream os(filename.c_str());

	copy(
		data,data+centersMap.num_elements(),
		ostream_iterator<float>(os, "\t")
		);
	os.close();
}

/**
    \brief Fourier transform, apply Mask, inverse FFT and normalize the result between 0 and 255.
    The heaviest part of the algorithm
*/
void Tracker::FFTapplyMask()
{
    //boost::progress_timer ptimer;
    if(!quiet) cout << "FFT ... ";
    fftwf_execute(forward_plan);
    if(view) displaySpectrum();

    if(!quiet) cout << "applying filter ... ";
    //the complex type given by fftw has no multiplication operator (it's pure C)
    //but it is equivalent to float[2], two consecutive floats in memory : (real, imaginary)
    float *c = data, *last = data + 2*FFTmask.num_elements();
    bool *b = FFTmask.origin();
    while(c!=last)
    {
    	*c++ *= *b;
    	*c++ *= *b++;
    }
    if(view) displaySpectrum("Filtered Spectrum");

    if(!quiet) cout << "inverse FFT ... ";
    fftwf_execute(backward_plan);
}

/** \brief Function object class for less-than inequality comparison between pointer*/
template <class T> struct less_pointer : binary_function <T*,T*,bool> {
  bool operator() (const T* x, const T* y) const
    {return *x < *y;}
};

template <class T> struct deref : unary_function<T*, T> {
    bool operator()(const T* x) const
        {return *x;} 
};

/** @brief Fill the centersMap so that only pixels corresponding to a center are True
	Works like a dilatation-origin and a threshold (see Peter Lu) but with no extra memory consumption.
	Hopefully this is also faster (~3N instead of N^3)
*/
void Tracker::findPixelCenters(float threshold)
{
    //boost::progress_timer ptimer;
    const size_t margin = 6;
	//We work only with unpadded, normalized images
	unpad();
	//normalize();
	//rather than normalizing the whole image, un-normalizing the threshold to be more efficient
	//std::pair<float*, float*> minmax = boost::minmax_element(data, data+centersMap.num_elements());
	//cout<< *minmax.first <<"<"<< *minmax.second<<endl;
	//threshold = threshold *(*minmax.second - *minmax.first)/255.0 + *minmax.first;
	threshold *= centersMap.num_elements();
	/*double sumI = accumulate(data, data+centersMap.num_elements(), 0.0, plus<double>());
	if(!hasRefBrightness)
	{
        refBrightness = sumI;
        hasRefBrightness = true;
	}
    else
        threshold += (sumI-refBrightness);*/

	//what is the offest of the first usefull pixel ?
	size_t offset =0;
	for(size_t d=0;d<3;++d)
        if(centersMap.shape()[d] >= 2*margin+1)
            offset += margin * centersMap.strides()[d];
    //pointer to the first useful pixel value
    float *px = data + offset, *end = data+centersMap.num_elements()-offset;

	//what are the neighbours of the first usefull pixel ? (pointers)
	vector<float*> neighbours;
	for(size_t d=0;d<3;++d)
        if(centersMap.shape()[d] >= 2*margin+1)
        {
            neighbours.push_back(px-centersMap.strides()[d]);
            neighbours.push_back(px+centersMap.strides()[d]);
        }
    
    //pointer to the center of not value of the first useful pixel
	bool *c = centersMap.origin() + offset;
	
	//what are the preceding potential neighbours of the first useful pixel ? (pointers)
    vector<bool*> ngbcenters;
    for(size_t d=0;d<3;++d)
        if(centersMap.shape()[d] >= 2*margin+1)
            ngbcenters.push_back(c-centersMap.strides()[d]);
	
	//iterate over the useful pixels, incrementing the neighbourhoods pointers
	while(px != end)
	{
        //set center to false if px<threshold || px< one of its neighbours
        //set center to true if px>=threshold && !(px< one of its neighbours)
        //cout<<"px="<< *px<<"="<<more_than_pointer<float>(*px).big<<endl;
        //exit(0);
        *c++ =
        (
            (*px >= threshold) &&
            (
                find_if(
                    neighbours.begin(), neighbours.end(),
                    bind1st(less_pointer<float>(), px)
                    ) == neighbours.end()
                ) &&
            (//Remove the adjacent neighbours (two pixels exactly equal)
                find_if(
                    ngbcenters.begin(), ngbcenters.end(),
                    deref<bool>()
                    ) == ngbcenters.end()
                )
        );
        //incrementation
        px++;
        for(vector<float*>::iterator n = neighbours.begin(); n!=neighbours.end(); ++n)
            (*n)++;
        for(vector<bool*>::iterator nc = ngbcenters.begin(); nc!=ngbcenters.end(); ++nc)
            (*nc)++;
	}

    //removing margin in each dimension
    /*if(centersMap.shape()[0] >= 2*margin+1)
    {
        fill(centersMap.origin(), centersMap.origin()+centersMap.strides()[0]*margin, false);
        fill(centersMap.origin()+centersMap.num_elements()-centersMap.strides()[0]*margin, centersMap.origin()+centersMap.strides()[0]*margin, false);
    }*/
    if(centersMap.shape()[1] >= 2*margin+1)
        for(boost::multi_array<bool,3>::iterator i=centersMap.begin(); i<centersMap.end();++i)
        {
            fill_n(i->origin(), i->strides()[0]*margin, false);
            fill_n(i->origin() + (i->size()-margin)*i->strides()[0], i->strides()[0]*margin, false);
        }
    if(centersMap.shape()[2] >= 2*margin+1)
        for(boost::multi_array<bool,3>::iterator i=centersMap.begin(); i<centersMap.end();++i)
            for(boost::multi_array<bool,3>::subarray<2>::type::iterator j=i->begin(); j<i->end();++j)
            {
                fill_n(j->begin(), margin, false);
                fill_n(j->rbegin(), margin, false);
            }
            
}

/** @brief compare two indicies by the intensity of the pixel they point to  */
struct compIntensities : public binary_function<const size_t&, const size_t&, float>
{
	float *data;
	explicit compIntensities(float *d) {data=d;};
	bool operator()(const size_t& i, const size_t& j)
	{
		return *(data + i) < *(data + j);
	};
};

/** \brief get the centroid of the neighbourhood of an image pixel given by it's offset */
struct centroid : public std::unary_function<const size_t&, std::valarray<double> >
{
    const boost::multi_array_ref<float,3> &image;
    boost::multi_array<float,3> ngb;
    //boost::multi_array<char, 4> slides;
    //first and second order derivative on 5 points
    static const boost::array<double,5> coefprime, coefsec;

    centroid(const boost::multi_array_ref<float,3> & im);//:image(im){return;};

    valarray<double> operator()(const size_t& l) const;
};
const boost::array<double,5> 
    centroid::coefprime  = {{1,-8, 0, 8, -1}}, 
    centroid::coefsec = {{-1, 16, -30, 16, -1}};

/** @brief centroid functor constructor. prepare the gradients  */
centroid::centroid(const boost::multi_array_ref<float,3> & im) : image(im)
{
    boost::array<size_t, 3> kernel_dims;
    for(size_t d=0;d<3;++d)
        kernel_dims[d] = (im.shape()[d]<3)?1:3;
    
    ngb.resize(kernel_dims);

    
}


/** \brief get the centroid of the neighbourhood of an image pixel given by it's offset */
valarray<double> centroid::operator()(const size_t& l) const
{
    const int scope = 2;
    //convert the raveled index to 3D indices
	size_t
		i = l / image.strides()[0],
		j = (l % image.strides()[0]) / image.strides()[1],
		k = (l % image.strides()[0]) % image.strides()[1];
    
    //the data of the neighbourhood view are copied together for the coder's sanity
	boost::multi_array<float,3> ngb =
		image[boost::indices
				[image.shape()[0]<2*scope+1 ? range() : range(i-scope, i+scope+1)]
				[image.shape()[1]<2*scope+1 ? range() : range(j-scope, j+scope+1)]
				[image.shape()[2]<2*scope+1 ? range() : range(k-scope, k+scope+1)]
			];
    //Find the extrema of the neighbourhood.
    std::pair<float*, float*> minmax = boost::minmax_element(ngb.origin(), ngb.origin()+ngb.num_elements());
    //If the neighbourhood contains a negative pixel, we are at the edge of a Fourier filtering artefact that should not be considered a particle
    if(*minmax.first < 0)
        return valarray<double>(-1.0, 3);
	//marking non local maxima (including diagonals)
	if(image.origin()[l] != *minmax.second)
        return valarray<double>(-1.0, 3);
    
    valarray<double> c(0.0,3);
    
    double d =0, dd=0;
    for(int mi=0; mi<2*scope+1; mi++)
    {
        d += (double)ngb[scope][mi][scope] * coefprime[mi];
        dd += (double)ngb[scope][mi][scope] * coefsec[mi];
        //std::cout<<  ngb[scope][mi][scope] << ", "; //TODO: export to an file
    }
    //std::cout<<std::endl;
    //exit(0);
    c[1] = - d/dd;
    //c[1] = 0.4383 - d/dd;
    
    
    
    
    c[0] += i;
	c[1] += j;
	c[2] += k;
	//c[0] = d;
	//c[2] = dd;
	return c;
    

	

	/*//calculation of the intensity centroid
	valarray<double> c(0.0,3);
	double total_w = 0.0;
	float *px = ngb.origin();
	for(int x=0; x<ngb.shape()[0];++x)
        for(int y=0; y<ngb.shape()[1];++y)
            for(int z=0; z<ngb.shape()[2];++z)
            {
                const double weight = pow((double)(x-scope), 2) + pow((double)(y-scope), 2) + pow((double)(z-scope), 2) * (double)(*px);
                c[0] += (x-scope)*weight;
                c[1] += (y-scope)*weight;
                c[2] += (z-scope)*weight;
                total_w += weight ;
                px++;
            }

    c /= total_w/pow(2.0*scope+1, 2);
	c[0] += i;
	c[1] += j;
	c[2] += k;
	return c;*/
};

/** \brief helper functor to transform zyx into xyz*/
template<class Coordinates>
struct revert : public unary_function<Coordinates& ,void>
{
	void operator()(Coordinates& v){swap(v[0], v[2]);}
};

template<class Coordinates>
struct markedPositionRemover : public unary_function<Coordinates& ,bool>
{
	bool operator()(Coordinates& v){return v[0]<0;}
};


/** @brief get the centers as intensity centroids of local maxima neighbourhood
  *
  * Local maxima are taken as the bright pixels of centersMap
  */
Particles Tracker::getSubPixel(bool autoThreshold)
{
    //boost::progress_timer ptimer;
	//get the positions of the bright centers from flatten centersMap
	if(!quiet) cout<<"get positions ... ";
	list<size_t> flat_positions;
	bool* b = centersMap.origin();
	for(size_t i=0; i<centersMap.num_elements(); ++i)
		if(*b++)
			flat_positions.push_back(i);

	//sort by increasing pixel intensity
	if(!quiet) cout<<flat_positions.size()<<" potential centers ... ";
	flat_positions.sort(compIntensities(data));

	//Find the best intensity threshold
	if(autoThreshold)
	{
	    //construct the centers intensity histogram
	    const long double minhist = *(data+flat_positions.front()),
            maxhist = *(data+flat_positions.back());
        const long double scale = 100.0/(maxhist-minhist);
	    vector<size_t> hist(100, 0);
	    for(list<size_t>::const_iterator c = flat_positions.begin(); c!= flat_positions.end(); ++c)
	    {
	    	size_t l = (*(data + (*c)) - minhist) * scale;
	    	if(l<hist.size())
	    		hist[l]++;
	    }


        //find the population global maximum
        vector<size_t>::iterator peak = max_element(hist.begin(), hist.end());
        if(!quiet) cout<<"peak ="<<distance(hist.begin(), peak)<<" ... ";

        //find the last non-null bin before the peak, and at the same time the bin with the minimum population between it and the peak
        vector<size_t>::iterator non_z = lower_bound(hist.begin(), peak, 1);
        vector<size_t>::iterator lowest = min_element(non_z, peak);
        //cout<<"non_z="<<distance(hist.begin(),non_z)<<"\tlowest="<<distance(hist.begin(),lowest)<<endl;
        while((*lowest)==0)
        {
            non_z = lowest;
            lowest = min_element(lowest, peak);
            //cout<<"non_z="<<distance(hist.begin(),non_z)<<"\tlowest="<<distance(hist.begin(),lowest)<<endl;
        }

        //lowest gives the position of the valley of the histogram, ie the best threshold
        if(!quiet) cout<<"best threshold = "<<distance(hist.begin(), lowest)/scale + minhist<<" ("<<distance(hist.begin(), lowest)<<"%) ... ";
        //We remove all centers that have lower intensity than the best threshold
        list<size_t>::iterator i = flat_positions.begin();
        advance(i, accumulate(hist.begin(), lowest, (size_t)0));
        flat_positions.erase(flat_positions.begin(), i);
        if(!quiet) cout<<flat_positions.size()<<" centers after auto thresholding ... ";
	}

	//centroid binning
	if(!quiet) cout<<"centroid ... ";
	list< valarray<double> > positions;
	boost::array<size_t,3> shape;
	copy(centersMap.shape(), centersMap.shape()+3, shape.begin());
    transform(
		flat_positions.rbegin(), flat_positions.rend(),
		back_inserter(positions),
		centroid(boost::multi_array_ref<float,3>(data,shape))
		);
    //removing marked centers
    positions.remove_if(markedPositionRemover<valarray<double> >());
    if(!quiet) cout<<positions.size()<<" are real local maxima ... ";

    Particles centers;
	for(size_t d=0;d<3;++d)
	{
        centers.bb.edges[d].first  = 0;
		centers.bb.edges[d].second = centersMap.shape()[d];
	}
	//centers.reserve(positions.size());
	copy(positions.begin(), positions.end(), back_inserter(centers));

	if(fortran_order)
	{
		//If the image was filled in row major ordering of the dimensions, the coordinates are given as z,y,x by the tracker.
		//Here we convert to the human-readable x,y,z coordinates
		for_each(centers.begin(), centers.end(), revert<valarray<double> >());
		//don't forget the bounding box
		swap(centers.bb.edges[0].second, centers.bb.edges[2].second);
	}
    if(!quiet) cout << "done!" << endl;
    if(view)
    {
        markCenters();
        display("Centers");
    }
	return centers;
}

struct centerMarker : public binary_function<const float&, const bool&, float>
{
    float v;
    centerMarker(float value) {v=value;}
    float operator()(const float& x, const bool& b){return (b?v:x);}
};

/** @brief Mark the centers on the image. Be carefull, data is overwritten. */
void Tracker::markCenters()
{
    std::pair<float*, float*> minmax = boost::minmax_element(data, data+centersMap.num_elements());
    transform(
        data, data+centersMap.num_elements(),
        centersMap.origin(), data,
        centerMarker((*minmax.second - *minmax.first))
        );
}

/** @brief Free the memory and leave the iterator as the default constructed one. Tracker get deleted  */
void TrackerIterator::close()
{
    if(tracker) delete tracker;
    if(centers) delete centers;
}

/** @brief set the band pass filter in order to select the same real sizes in x, y and z directions  */
void TrackerIterator::setIsotropicBandPass(double radiusMin, double radiusMax)
{
    boost::array<double,3>
        radiiMin = {{radiusMin/getZXratio(), radiusMin, radiusMin}},
        radiiMax = {{radiusMax/getZXratio(), radiusMax, radiusMax}};
    this->tracker->makeBandPassMask(radiiMin, radiiMax);
}

/** @brief custom setting of the band pass filter. Use at your own risk.  */
void TrackerIterator::setAnisotropicBandPass(double radiusMin, double radiusMax, double zRadiusMin, double zRadiusMax)
{
    boost::array<double,3>
        radiiMin = {{zRadiusMin, radiusMin, radiusMin}},
        radiiMax = {{zRadiusMax, radiusMax, radiusMax}};
    this->tracker->makeBandPassMask(radiiMin, radiiMax);
}

/** @brief set the band pass filter in order to select the same real sizes in x, y and z directions  */
void TrackerIterator::setIsotropicLowPass(double radiusMin)
{
    boost::array<double,3>
        radiiMin = {{radiusMin/getZXratio(), radiusMin, radiusMin}};
    this->tracker->makeLowPassMask(radiiMin);
}

/** @brief get the file name to export the intensity to at the present time step.  */
string TrackerIterator::getIntensityFile()
{
    assert(!!this->IntensitySerie);
    return (*this->IntensitySerie)%this->time_step;
}




/** @brief track particles and return the result.
  *
  * Contrary to the dereferencement operator for standard iterator, this is a heavy operation.
  * The result is cached until the next incrementation or setTimeStep
  */
Particles& TrackerIterator::operator*()
{
    if(!centers)
    {
        if(reachedEnd())
        {
            cout<<"reached end"<<endl;
            centers = new Particles();
        }
        else
        {
            //Get the coordinate of the centers expressed in pixel units
            if(this->noThreshold)
                centers = new Particles(tracker->trackXYZ(tracker->mean, true));
            else
                centers = new Particles(tracker->trackXYZ(this->threshold));


            //export intensities if asked
            if(!!this->IntensitySerie)
            {
                if(!getTracker().quiet) cout <<"export intensities"<<endl;
                ofstream out(getIntensityFile().c_str());
                vector<float> intensities = tracker->getIntensities(*centers);
                copy(
                    intensities.begin(), intensities.end(),
                    ostream_iterator<float>(out,"\n")
                    );
            }

            //The real size of the pixel in z is not in general the same as the size in x or y
            //conversion to real size
            valarray<double> v(1.0,3);
            v[2] = getZXratio();
            (*centers) *= v;
            if(!getTracker().quiet) cout <<"z multiplied by "<<v[2]<<endl;
        }
    }
    return *centers;
}

