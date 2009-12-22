#include <CImg.h>
using namespace cimg_library;

int main(int argc, char* argv[]) {
  const CImg<unsigned char> img = CImg<unsigned char>(argv[1]).channel(1);
  CImgList<> F = img.get_FFT();
  cimglist_apply(F,translate)(img.dimx()/2,img.dimy()/2,0,0,2);
  const CImg<unsigned char> mag = ((F[0].get_pow(2) + F[1].get_pow(2)).sqrt()+1.0f).log().normalize(0,255);
  CImgList<unsigned char> visu(img,mag);
  CImgDisplay disp(visu,"[#16] - Fourier Filtering (Click to set filter)");
  CImg<unsigned char> mask(img.dimx(),img.dimy(),1,1,1);
  unsigned char one[] = { 1 }, zero[] = { 0 }, white[] = { 255 };
  int rmin = 0, rmax = 256;
  while (!disp.is_closed && !disp.is_keyQ && !disp.is_keyESC) {
    disp.wait();
    const int
      xm = disp.mouse_x*2*img.dimx()/disp.dimx()-img.dimx(),
      ym = disp.mouse_y*img.dimy()/disp.dimy(),
      x = xm-img.dimx()/2,
      y = ym-img.dimy()/2;
    if (disp.button && xm>=0 && ym>=0) {
      const int r = (int)cimg::max(0.0f,(float)std::sqrt((float)x*x+y*y)-3.0f);
      if (disp.button&1) rmax = r;
      if (disp.button&2) rmin = r;
      if (rmin>=rmax) rmin = cimg::max(rmax-1,0);
      mask.fill(0).draw_circle(mag.dimx()/2,mag.dimy()/2,rmax,one).
        draw_circle(mag.dimx()/2,mag.dimy()/2,rmin,zero);
      CImgList<> nF(F);
      cimglist_for(F,l) nF[l].mul(mask).translate(-img.dimx()/2,-img.dimy()/2,0,0,2);
      visu[0] = nF.FFT(true)[0].normalize(0,255);
    }
    if (disp.is_resized) disp.resize(disp.window_dimx(),disp.window_dimx()/2).display(visu);
    visu[1] = mag.get_mul(mask).draw_text(5,5,white,zero,11,0.6f,"Freq Min/Max = %d / %d",(int)rmin,(int)rmax);
    visu.display(disp);
  }
  return 0;
}
