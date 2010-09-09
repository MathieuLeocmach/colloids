/**
    Copyright 2010 Mathieu Leocmach

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

 * \file radiiTracker.hpp
 * \brief Define classes to extract particles radii from 3D image data
 * \author Mathieu Leocmach
 *
 *
 */


#ifndef radii_tracker_H
#define radii_tracker_H

extern "C" {
#include "fftw3.h"
}

#include "particles.hpp"
#include "lifFile.hpp"
#include "files_series.hpp"
#include <boost/multi_array.hpp>
namespace Colloids{

/** \brief basic radii tracker class containing the radii extraction algorythm*/
class RadiiTracker
{
    public:
        typedef unsigned char                                       pixel_type;
		typedef boost::multi_array<pixel_type, 3>					array_type;

		//size_t channel;
		bool view, quiet, fortran_order;
		double ZXratio, rmin;

        explicit RadiiTracker(const boost::array<size_t,3> &dims, const unsigned fs=FFTW_ESTIMATE) : flags(fs), data(dims)
        {
            ZXratio = 1.0;
            precision = 0.1;
            rmin = 1.0;
            rmax = 10.0;
            resize();
        };
		~RadiiTracker()
		{
		    fftwf_destroy_plan(forward);
            fftwf_destroy_plan(backward);
            fftwf_free(in); fftwf_free(out);
        }

        template <class InputIterator>
        InputIterator fillImage_charToUchar(InputIterator first);
        template <class OutputIterator>
        OutputIterator copyImage(OutputIterator result) const;

        Coord getNeighbourhood(const Coord &point, array_type &ngb) const;
        double getRadius(const array_type &ngb, const Coord &center, const double &rmax);
        std::vector<double> getRadii();

        void setCenters(Particles c);
        void setPrecision(const double &p, bool noresize=false){this->precision=p; resize();}
        void setRmax(const double &r, bool noresize=false){this->rmax=r; resize();}
        double getPrecision(){return precision;}
        double getRmax(){return rmax;}

    protected:
		unsigned flags;
		array_type data;
		Particles centers;
		double precision, rmax;
		int N;
		fftwf_plan forward, backward;
		float *in;
        fftwf_complex *out;

        void resize();
};

/** \brief Virtual glue between RadiiTracker and data source
Can be used as a stadard InputIterator
*/
class RadiiTrackerIterator : public std::iterator<std::input_iterator_tag, std::vector<double> >
{
 	protected:
        RadiiTracker *tracker;
        std::vector<double> *radii;
        size_t channel, time_step;
        FileSerie *centersSerie;

 	public:
		/** \brief default constructor that should be used only to get an "end" iterator*/
		RadiiTrackerIterator(){};
		virtual ~RadiiTrackerIterator(){};
		virtual void close();

		const RadiiTracker& getTracker() const {return *tracker;};
		void setCentersSerie(const FileSerie &cSerie){centersSerie = new FileSerie(cSerie);}
		void setCenters(Particles c);

		void setView(bool v){tracker->view=v;}
		bool view() const {return tracker->view;}
		void setQuiet(bool v){tracker->quiet=v;}
		bool quiet() const {return tracker->quiet;}

		void setChannel(size_t ch){channel=ch;};
		size_t getChannel(){return channel;};
		virtual void setTimeStep(size_t t)=0;
		size_t getTimeStep(){return this->time_step;};
		virtual double getZXratio()=0;

		void setPrecision(const double &p){tracker->setPrecision(p);}
		double precision() const {return tracker->getPrecision();}
		void setRmin(const double &r){tracker->rmin=r;}
		double rmin() const {return tracker->rmin;}
		void setRmax(const double &r){tracker->setRmax(r);}
		double rmax() const {return tracker->getRmax();}

		virtual RadiiTrackerIterator& operator++()=0;
		virtual std::vector<double>& operator*();
		bool operator!=(const RadiiTrackerIterator& rhs) {return time_step!=rhs.time_step;}
};

/** \brief glue between RadiiTracker and LifFile
Can be used as a stadard InputIterator
*/
class LifRadiiTracker : public RadiiTrackerIterator
{
 	LifSerie *serie;
 	std::istreambuf_iterator<char> iterator;

 	public:
		explicit LifRadiiTracker(LifSerie &serie, const std::string &centerSeriePrefix, const size_t ch=0, const unsigned fs=FFTW_ESTIMATE);
		/** \brief default constructor that should be used only to get an "end" iterator*/
		LifRadiiTracker end(){return LifRadiiTracker(getLif().getNbTimeSteps());};

		LifSerie& getLif() const {return *serie;};

		void chooseChannel();
		void setTimeStep(size_t t);
		double getZXratio(){return serie->getZXratio();};

		LifRadiiTracker& operator++();
		bool operator==(const LifRadiiTracker& rhs) {return time_step==rhs.time_step && serie==rhs.serie;}
		//bool operator!=(const LifTracker& rhs) {return time_step!=rhs.time_step;}


	private:
        LifRadiiTracker(const size_t end_step) : iterator(){time_step = end_step;};
		boost::array<size_t,3> getTrackerDims() const;
		std::streampos tellg(){return serie->tellg();}

};

/** @brief fill the Image from iterators thinking they are contening char when their content is in fact unsigned char
	Reading from a raw unsigned char image
	\code
	ifstream in("C:/Code_data/img.raw", ios_base::in | ios_base::binary);
	track2Dimg.fillImage_charToUchar(ordering,istreambuf_iterator<char>(in));
	in.close();
	\endcode
 */
template <class InputIterator>
InputIterator RadiiTracker::fillImage_charToUchar(InputIterator first)
{
	pixel_type *d = data.origin();
	for(size_t l=0; l<data.num_elements(); ++l)
        *d++ = static_cast<pixel_type>(*first++);
	return first;
}

/** @brief outputImage
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
OutputIterator RadiiTracker::copyImage(OutputIterator result) const
{
		return copy(data.origin(), data.origin()+data.num_elements(), result);
}
}
#endif
