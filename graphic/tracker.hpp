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

 * \file tracker.hpp
 * \brief Define classes to track particles in 3D image series
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 7 april 2009
 *
 * Object oriented implementation of the code elaborated in 2008.
 * Based on Crocker & Grier algorithm
 *
 */


#ifndef tracker_H
#define tracker_H

extern "C" {
#include "fftw3.h"
}
//#include "lifFile.hpp"
//#include "../files_series.hpp"
#include "particles.hpp"
//#include <CImg.h>
#include <boost/multi_array.hpp>
namespace Colloids{
/** \brief basic tracker class containing the tracking algorythm*/
class Tracker
{
    public:
		typedef boost::multi_array_ref<float, 3>					array_type_r;
		typedef boost::multi_array<fftwf_complex, 3>				array_type_c;
		typedef boost::multi_array<bool, 3>							array_type_b;
		typedef boost::multi_array<float,3>::array_view<3>::type	view_type;

		//size_t channel;
		bool view, quiet, fortran_order;
		array_type_b FFTmask, centersMap;
		//double displayRadius;
		//double Zratio;

		explicit Tracker(const boost::array<size_t,3> &dims, const unsigned fs=FFTW_ESTIMATE);
		~Tracker();

		void setFlags(const unsigned fs=FFTW_ESTIMATE){flags = fs;};
        void setDimensions(const boost::array<size_t,3> &dims);
        void makeBandPassMask(const boost::array<double,3> &radiiMin, const boost::array<double,3> &radiiMax);
        void FFTapplyMask();
        void findPixelCenters(float threshold=0);

        array_type_r get_padded_image(const boost::array<size_t,3>& ordering);
        view_type get_image(const boost::array<size_t,3>& ordering );

        void unpad();
        void normalize();
        /*void setRefBrightness();
        void unsetRefBrightness(){hasRefBrightness=false;};*/

        Particles trackXYZ(const float &threshold=0.0f);
        std::vector<float> getIntensities(const Particles &centers);

        //std::vector<Particles*> granulometry(const double &radiusMin, const double &radiusMax);
        static bool loadWisdom(const std::string & filename);
		static void saveWisdom(const std::string & filename);
		template <class InputIterator>
        InputIterator fillImage(InputIterator first);
        template <class InputIterator>
        InputIterator fillImage_charToUchar(InputIterator first);
        template <class InputIterator>
        InputIterator fillSlice(const size_t slice, InputIterator first);
        template <class OutputIterator>
        OutputIterator copyImage(OutputIterator result) const;
        template <class OutputIterator>
        OutputIterator copyPaddedImage(OutputIterator result) const;
        template <class OutputIterator>
        OutputIterator copySpectrum(OutputIterator result) const;
		void maskToFile(std::string &filename) const;
		void imageToFile(std::string &filename) const;

		/** Display routines */
		void display(const std::string &windowName="Image") const;
		void displaySlice(const size_t dim, const size_t position, const std::string &windowName="") const;
		void displaySpectrum(const std::string &windowName="Spectrum") const;
		void displayMask(const std::string &windowName="Mask") const;
		void displayCenters(const std::string &windowName="Centers") const;
		void markCenters();

	private:
		unsigned flags;
		//array_type_b FFTmask, centersMap;
		float *data;
		fftwf_plan forward_plan,backward_plan;
		bool padded;/*, hasRefBrightness;
		double refBrightness;*/

		//void FFTapplyMask();
		//void findPixelCenters(const float &threshold=0);
		Particles getSubPixel();


		//bool compIntensities(size_t i, size_t j);


};

/** \brief Virtual glue between Tracker and data source
Can be used as a stadard InputIterator
*/
class TrackerIterator : public std::iterator<std::input_iterator_tag, Particles>
{
 	protected:
        Tracker *tracker;
        Particles *centers;
        size_t channel, time_step;
        float threshold;
        FileSerie *IntensitySerie;

 	public:
		/** \brief default constructor that should be used only to get an "end" iterator*/
		TrackerIterator(){};
		virtual ~TrackerIterator(){};
		virtual void close();

		const Tracker& getTracker() const {return *tracker;};

		virtual void setIsotropicBandPass(double radiusMin, double radiusMax);
		virtual void setAnisotropicBandPass(double radiusMin, double radiusMax, double zRadiusMin, double zRadiusMax);
		void setThreshold(const float thr=0.0f) {this->threshold=thr;};
		const float& getThreshold() const {return this->threshold;};
		void setIntensitySerie(FileSerie &is){this->IntensitySerie=&is;};
		std::string getIntensityFile();
		void setView(bool v){tracker->view=v;}
		bool view() const {return tracker->view;}
		void setQuiet(bool v){tracker->quiet=v;}
		bool quiet() const {return tracker->quiet;}

		void setChannel(size_t ch){channel=ch;};
		size_t getChannel(){return channel;};
		virtual void setTimeStep(size_t t)=0;
		size_t getTimeStep(){return this->time_step;};
		virtual bool reachedEnd()=0;
		virtual double getZXratio()=0;

		virtual TrackerIterator& operator++()=0;
		virtual Particles& operator*();
		bool operator!=(const TrackerIterator& rhs) {return time_step!=rhs.time_step;}
};

/** @brief fillImage (add the necessaray padding)
	Reading from formatted ascii input
	\code
	ifstream in("C:/Code_data/img.txt");
	track2Dimg.fillImage(ordering,istream_iterator<float>(in));
	in.close();
	\endcode
*/
template <class InputIterator>
InputIterator Tracker::fillImage(InputIterator first)
{
	float *d = data;
	for(size_t i=0; i<centersMap.shape()[0];++i)
		for(size_t j=0; j<centersMap.shape()[1];++j)
		{
			for(size_t k=0; k<centersMap.shape()[2];++k)
				*d++ = *first++;
			//padding the last dimension
			for(size_t k=centersMap.shape()[2]; k< 2 * FFTmask.shape()[2];++k)
				d++;
		}
	this->padded = true;
	return first;
}

/** @brief fill the Image from iterators thinking they are contening char when their content is in fact unsigned char
	Reading from a raw unsigned char image
	(add the necessaray padding)
	\code
	ifstream in("C:/Code_data/img.raw", ios_base::in | ios_base::binary);
	track2Dimg.fillImage_charToUchar(ordering,istreambuf_iterator<char>(in));
	in.close();
	\endcode
 */
template <class InputIterator>
InputIterator Tracker::fillImage_charToUchar(InputIterator first)
{
    //using namespace boost::accumulators;
    //accumulator_set<long double, features<tag::min, tag::max, tag::mean> > acc;
	float *d = data;
	for(size_t i=0; i<centersMap.shape()[0];++i)
		for(size_t j=0; j<centersMap.shape()[1];++j)
		{
			for(size_t k=0; k<centersMap.shape()[2];++k)
				*d++ = static_cast<unsigned char>(*first++);

			//padding the last dimension
			for(size_t k=centersMap.shape()[2]; k< 2 * FFTmask.shape()[2];++k)
				*d++ = 0.0f;
		}
	this->padded = true;
	//std::cout<< " original image min="<<min(acc)<<" max="<<max(acc)<<" sum="<<mean(acc)<<" ... first pixel: "<<*data<<" ... ";
	return first;
}

/** @brief fill a slice of the image (add the necessaray padding)*/
template <class InputIterator>
InputIterator Tracker::fillSlice(const size_t slice, InputIterator first)
{
	float *d = data + slice*FFTmask.strides()[0]*2;
    for(size_t j=0; j<centersMap.shape()[1];++j)
    {
        for(size_t k=0; k<centersMap.shape()[2];++k)
            *d++ = *first++;
        //padding the last dimension
        for(size_t k=centersMap.shape()[2]; k< 2 * FFTmask.shape()[2];++k)
            d++;
    }
	this->padded = true;
	return first;
}

/** @brief outputImage (removing the paddings)
	Writting to formatted ascii
	\code
		ofstream out("C:/Code_output/img.txt", ios_base::out | ios_base::trunc);
		track2Dimg.copyImage(ostream_iterator<float>(out,"\t"));
		out.close();
	\endcode
	Writting to binary
	\code
		ofstream out("C:/Code_output/img.raw", ios_base::out | ios_base::trunc | ios_base::binary);
		track2Dimg.copyImage(ostreambuf_iterator<char>(out));
		out.close();
	\endcode
*/
template <class OutputIterator>
OutputIterator Tracker::copyImage(OutputIterator result) const
{
	if(this->padded)
	{
		float *d = data;
		for(size_t i=0; i<centersMap.shape()[0];++i)
			for(size_t j=0; j<centersMap.shape()[1];++j)
			{
				for(size_t k=0; k<centersMap.shape()[2];++k)
					*result++ = *d++;
				//padding the last dimension
				for(size_t k=centersMap.shape()[2]; k<2*FFTmask.shape()[2];++k)
					*d++ = 0.0f;
			}
		return result;
	}
	else
		return copy(data, data+centersMap.num_elements(), result);
}

/** @brief output Image (with the paddings)
*/
template <class OutputIterator>
OutputIterator Tracker::copyPaddedImage(OutputIterator result) const
{
	return copy(data, data+2*FFTmask.num_elements(), result);
}

/** @brief output Spectrum
*/
template <class OutputIterator>
OutputIterator Tracker::copySpectrum(OutputIterator result) const
{
	float *d = data;
	for(size_t l=0; l<FFTmask.num_elements();++l)
		*result++ = log(1+sqrt(pow(*d++, 2.0) + pow(*d++, 2.0)));
	return result;
}
}
#endif
