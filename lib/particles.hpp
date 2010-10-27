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

 * \file particles.hpp
 * \brief Defines classes for particles
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 26 Novemeber 2008
 *
 * Define the mother of all Particles classes
 *
 */


#ifndef particles_H
#define particles_H

#include "index.hpp"
#include "fields.hpp"
#include "boo_data.hpp"

#include <boost/multi_array.hpp>
#include <boost/bind.hpp>

#include <algorithm>
#include <deque>
#include <valarray>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include <stdexcept>

//#include <tvmet/Vector.h>



namespace Colloids
{
    typedef RStarIndex_S::RTree                     RTree;
    typedef std::vector< std::vector<size_t> >         NgbList;

    struct Bond : private std::pair<size_t, size_t>
	{
		explicit Bond(const size_t &x, const size_t &y){this->assign(x,y);}
		Bond(std::pair<size_t, size_t> p){this->assign(p.first, p.second);};
		Bond(){first=static_cast<size_t>(-2); second=static_cast<size_t>(-1);}

		void assign(const size_t &x, const size_t &y)
		{
			if(x<y)
			{
				this->first=x;
				this->second=y;
			}
			else
			{
				this->first=y;
				this->second=x;
			}
		};
		const size_t& low() const {return this->first;}
		const size_t& high() const {return this->second;}
		bool operator<(const Bond &rhs) const
		{
			return (this->first < rhs.first) || (this->first == rhs.first && this->second < rhs.second);
		}
	};

    typedef std::set<Bond>	BondSet;

    BondSet ngb2bonds(const NgbList& ngbList);

    /**
        \brief defines a set of particles having the same radius
    */
    class Particles : public std::vector<Coord>
    {
        /** \brief A spatial index of the particles */
        std::auto_ptr<SpatialIndex> index;

        /** \brief A neighbour list */
        std::auto_ptr<NgbList> neighboursList;

        public:
            /** \brief overall bounding box */
            BoundingBox bb;
            /** \brief (mean) radius of all the particles */
            double radius;


            /** \brief constructors and destructor */
            Particles(void) : std::vector<Coord>(0,Coord(0.0,3)){radius=1.0;return;};
            Particles(const std::vector<Coord> &data, const double &r=1.0) : std::vector<Coord>(data){radius=r;};
            Particles(const size_t &n, const double &d=0.0, const double &r=1.0);
            Particles(const std::string &filename, const double &r=1.0);
            Particles(const size_t &Nb, const BoundingBox &b, const std::string &filename, const double &r=1.0);
            virtual ~Particles(){return;}

            void push_back(const Coord &p);

            Particles cut(const double &sep);

            /** Geometric transforms    */
            Particles& operator*=(const Coord &v);
            Particles& operator*=(const double &mul);
            virtual Particles& operator+=(const Coord &v);

            /** Geometry related */
            virtual Coord getDiff(const Coord &from,const size_t &to) const;
            virtual Coord getDiff(const size_t &from,const size_t &to) const;
            virtual double getAngle(const size_t &origin,const size_t &a,const size_t &b) const;
            virtual std::vector<size_t> selectInside_noindex(const double &margin) const;
            void loadInside(std::vector<size_t> &inside) const;

            /** Index related   */
            static BoundingBox bounds(const Coord &center,const double &r=0.0);
            bool hasIndex() const {return !!index.get();};
            void setIndex(SpatialIndex *I){index.reset(I);};
            void makeRTreeIndex();
            BoundingBox getOverallBox() const;

            /** Spatial query and neighbours. Depends on both geometry and spatial index */
            virtual std::vector<size_t> selectEnclosed(const BoundingBox &b) const;
            std::vector<size_t> getEuclidianNeighbours(const Coord &center, const double &range) const;
            std::vector<size_t> getEuclidianNeighbours(const size_t &center, const double &range) const;
            size_t getNearestNeighbour(const Coord &center, const double &range=1.0) const;
            std::multimap<double,size_t> getEuclidianNeighboursBySqDist(const Coord &center, const double &range) const;
            NgbList & makeNgbList(const double &bondLength);
            NgbList & makeNgbList(const BondSet &bonds);
            const NgbList & getNgbList() const {return *this->neighboursList;};
            void delNgbList(){neighboursList.reset();};
            BondSet getBonds() const {return ngb2bonds(getNgbList());};
            virtual std::vector<size_t> selectInside(const double &margin) const;


            /**Bond Orientational Order related */
            BooData sphHarm_OneBond(const size_t &center, const size_t &neighbour) const;
            BooData getBOO(const size_t &center) const;
            BooData getCgBOO(const std::vector<BooData> &BOO, const size_t &center) const;
            void getBOOs(std::vector<BooData> &BOO) const;
            void getBOOs(const std::vector<size_t> &selection, std::vector<BooData> &BOO) const;
            void getCgBOOs(const std::vector<size_t> &selection, const std::vector<BooData> &BOO, std::vector<BooData> &cgBOO) const;
            void getSurfBOOs(std::vector<BooData> &BOO) const;
            void getBOOs_SurfBOOs(std::vector<BooData> &BOO, std::vector<BooData> &surfBOO) const;
            void getFlipBOOs(const std::vector<BooData> &BOO, std::vector<BooData> &flipBOO, const BondSet &bonds) const;
            void exportQlm(const std::vector<BooData> &BOO, const std::string &outputPath) const;
            void exportQ6m(const std::vector<BooData> &BOO, const std::string &outputPath) const;
            void load_q6m(const std::string &filename, std::vector<BooData> &BOO) const;
            void load_qlm(const std::string &filename, std::vector<BooData> &BOO) const;
            template<typename T> void removeOutside(const std::vector<size_t> &inside, std::vector<T> &BOO) const;

            /**Bond angle distribution related  */
            boost::array<double,180> getAngularDistribution(const size_t &numPt) const;
            boost::array<double,180> getMeanAngularDistribution(const NgbList &selection) const;

            /**Common neighbour analysis */
            bool is_ring(std::list<size_t> common) const;
            void getSP5c(std::vector< std::vector<size_t> > &SP5c) const;
            BondSet get1551pairs() const;
            BondSet get2331pairs() const;
            BondSet getSecondShell() const;

            /** histograms*/
            struct Binner : public std::binary_function<const size_t &,const size_t &,void>
            {
                const Particles & parts;
                size_t count;
                double cutoff;

                Binner(const Particles &p, const double &nbDiameterCutOff) : parts(p)
                {
                    count = 0;
                    cutoff = 2.0 * parts.radius * nbDiameterCutOff;
                };
                virtual ~Binner(void);
                virtual void operator()(const size_t &p, const size_t &q){};
                void operator<<(const std::vector<size_t> &selection);
            };

            struct RdfBinner : public Binner
            {
                std::vector<double> g;
                double scale;

                RdfBinner(const Particles &p, size_t n, const double &nbDiameterCutOff) : Binner(p,nbDiameterCutOff)
                {
                    g = std::vector<double>(n,0.0);
                    scale = n / cutoff;
                };
                void operator()(const size_t &p, const size_t &q)
                {
					g[(size_t)(norm2(parts.getDiff(p,q)) * scale)]++;
					count++;
				};
                void normalize(const size_t &n);
            };

            std::vector<double> getRdf(const std::vector<size_t> &selection, const size_t &n, const double &nbDiameterCutOff) const;
            std::vector<double> getRdf(const size_t &n, const double &nbDiameterCutOff) const;

            struct GlBinner : public RdfBinner
            {
                std::vector<double> gl;
                const std::vector<BooData> &boo;
                const size_t l;

                GlBinner(const Particles &p, size_t n, const double &nbDiameterCutOff, const std::vector<BooData> &BOO, const size_t &l)
                : RdfBinner(p,n,nbDiameterCutOff),boo(BOO), l(l)
                {
                    gl = std::vector<double>(n,0.0);
                };
                void operator()(const size_t &p, const size_t &q);
                void normalize(const size_t &n);
            };

            /** file outputs */
            void exportToFile(const std::string &filename) const;
            std::ostream & toVTKstream(std::ostream &out, const std::string &dataName = "particles") const;
            void exportToVTK(
                const std::string &filename,const BondSet &bonds,
                const std::vector<ScalarField> &scalars,	const std::vector<VectorField> &vectors,
                const std::string &dataName = "particles"
            ) const;
            void exportToVTK(const std::string &filename,
                const std::vector<ScalarField> &scalars,	const std::vector<VectorField> &vectors,
                const std::string &dataName = "particles"
            ) const;
            void exportToVTK(const std::string &filename,
                const std::vector<ScalarField> &scalars,
                const std::string &dataName = "particles"
            ) const;

            double getMinDim() const;
            virtual double getNumberDensity() const;
            double getVF() const;

            void loadBoo(const std::string &filename, boost::multi_array<double,2>&qw) const;
            //static bool areTooClose(const std::valarray<double> &c, const Coord &d,const double &Sep);

    };
    BondSet loadBonds(const std::string &filename);
    std::ostream &toVTKstream(std::ostream &out, const BondSet &bonds);
    inline std::ostream & operator<<(std::ostream& out, const Bond& b)
    {
    	out<<b.low()<<" "<<b.high();
    	return out;
    }
    inline std::istream & operator>>(std::istream& in, Bond& b)
    {
    	size_t x,y;
    	in>>x>>y;
    	b = Bond(x,y);
    	return in;
    }

    /**Inline functions, for performance*/

    /** \brief get the difference vector between a position and one of the particles */
    inline Coord Particles::getDiff(const Coord &from,const size_t &to) const
    {
        Coord diff(3);
        diff = at(to)-from;
        return diff;
    }

    /** \brief get the difference vector between two particles */
    inline Coord Particles::getDiff(const size_t &from,const size_t &to) const
    {
        Coord diff(3);
        diff = at(to)-at(from);
        return diff;
    }

    /** @brief get the indices of the particles enclosed by a query box  */
    inline std::vector<size_t> Particles::selectEnclosed(const BoundingBox &b) const
    {
        #ifndef NDEBUG
        if(!this->hasIndex()) throw std::logic_error("Set a spatial index before doing spatial queries !");
        #endif
        return (*index)(b);
    }

    /** @brief get the indices of the particles inside a reduction of the maximum bounding box  */
    inline std::vector<size_t> Particles::selectInside(const double &margin) const
    {
        #ifndef NDEBUG
        if(!this->hasIndex()) throw std::logic_error("Set a spatial index before doing spatial queries !");
        #endif
        return this->index->getInside(margin);
    }

    /**	\brief Bin a couple of particles into the g and g6 histogram. */
	inline void Particles::GlBinner::operator()(const size_t &p, const size_t &q)
	{
		if(!boo[p].isnull() && !boo[q].isnull())
		{
			count++;
			const size_t r = (size_t)(norm2(parts.getDiff(p, q)) * scale);
			g[r]++;
			gl[r] += boo[p].innerProduct(boo[q], l);
		}
	};

	/** @brief remove the values that are not in the selection      */
	template<typename T>
    void Particles::removeOutside(const std::vector<size_t> &inside, std::vector<T> &BOO) const
    {
        size_t p=0;
        for(std::vector<size_t>::const_iterator it = inside.begin(); it!=inside.end(); ++it)
        {
            while(p<*it)
                BOO[p++] = T();
            p=(*it)+1;
        }
        while(p<size())
            BOO[p++] = T();
    }



};
#endif

