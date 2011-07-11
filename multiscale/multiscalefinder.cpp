/*
 * multiscalefinder.cpp
 *
 *  Created on: 11 juil. 2011
 *      Author: mathieu
 */

#include "multiscalefinder.hpp"

using namespace std;

namespace Colloids {

MultiscaleFinder::MultiscaleFinder(const int nrows, const int ncols, const int nbLayers, const double &preblur_radius) :
		octaves(1, OctaveFinder(2*nrows, 2*ncols, nbLayers, preblur_radius))
{
	int ocr = nrows, occ = ncols;
	while(ocr >= 12 && occ >= 12)
	{
		this->octaves.push_back(OctaveFinder(ocr, occ, nbLayers, preblur_radius));
		ocr /= 2;
		occ /= 2;
	}

}

MultiscaleFinder::~MultiscaleFinder() {
	// TODO Auto-generated destructor stub
}

}
