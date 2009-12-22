#include <CImg.h>
#include <iostream>
using namespace cimg_library;

int main(int argc, char* argv[]) {
  const CImg<unsigned char> img = CImg<unsigned char>(argv[1]).channel(1).normalize(0,255);
  CImgList<> F = img.get_FFT();
  cimglist_apply(F,translate)(img.dimx()/2,img.dimy()/2,0,0,2);
  const CImg<unsigned char> mag = ((F[0].get_pow(2) + F[1].get_pow(2)).sqrt()+1.0f).log().normalize(0,255);
  CImg<unsigned char> tracked(img.dimx(),img.dimy(),img.dimz(),2,0);
  tracked.get_shared_channel(0)=img;
  CImg<unsigned char> treated(tracked);
  CImgList<float> visu(treated,mag,tracked);
  CImgDisplay disp(visu,"Fourier Filtering (Click on the middle image to set filter)");
  CImg<unsigned char> hist(256,65,1,1,0);
  CImgDisplay histo(hist, "Histogram (Clic to set threshold)");
  CImg<unsigned char> mask(img.dimx(),img.dimy(),1,1,1),thr_mask(mask);
  unsigned char one[] = { 1 }, zero[] = { 0 }, white[] = { 255 },grey[]={127};
  int fminx = 0, fmaxx = img.dimx(),fminy = 0, fmaxy = img.dimy();
  unsigned char thr=0;
  hist.draw_graph(img.get_histogram(),white).display(histo);
  while (!disp.is_closed && !disp.is_keyQ && !disp.is_keyESC) {
    disp.wait_all();
    const int
      xm = disp.mouse_x*3*img.dimx()/disp.dimx()-img.dimx(),
      ym = disp.mouse_y*img.dimy()/disp.dimy(),
      x = xm-img.dimx()/2,
      y = ym-img.dimy()/2;
    //histo.wait();
    const int xh = histo.mouse_x*hist.dimx()/histo.dimx();
    //std::cout << xh <<std::endl;
    if(histo.button && xh>=0 && xh<256)
    {
        thr=xh;
    }
    if (disp.button && xm>=0 && ym>=0 && xm<img.dimx() && ym < img.dimy()) {
      if(disp.is_keySHIFTLEFT || disp.is_keySHIFTRIGHT)
      {
            const int f = (int)cimg::max(0.0f,(float)std::sqrt((float)x*x+y*y));
            if (disp.button&1) fmaxx = fmaxy = f;
            if (disp.button&2) fminx = fminy = f;
      }
      else
      {
            if (disp.button&1)
            {
                if (!disp.is_keyCTRLRIGHT) fmaxx = std::abs(x);
                if (!disp.is_keyCTRLLEFT) fmaxy = std::abs(y);
            }
            if (disp.button&2)
            {
                if (!disp.is_keyCTRLRIGHT) fminx = std::abs(x);
                if (!disp.is_keyCTRLLEFT) fminy = std::abs(y);
            }
      }
      if (fminx>=fmaxx) fminx = cimg::max(fmaxx-1,0);
      if (fminy>=fmaxy) fminy = cimg::max(fmaxy-1,0);
      mask.fill(0).draw_ellipse(mag.dimx()/2,mag.dimy()/2,(float)fmaxx,(float)fmaxy,0,0,one).
        draw_ellipse(mag.dimx()/2,mag.dimy()/2,(float)fminx,(float)fminy,0,0,zero);
      CImgList<> nF(F);
      cimglist_for(F,l) nF[l].mul(mask).translate(-img.dimx()/2,-img.dimy()/2,0,0,2);
      visu[0].get_shared_channel(0) = nF.FFT(true)[0].normalize(0,255);
    }
    if (disp.is_resized) disp.resize(disp.window_dimx(),disp.window_dimx()/2).display(visu);
    visu[1] = mag.get_mul(mask).draw_text(5,5,white,zero,11,0.6f,"Radiusx Min/Max = %f / %f\nRadiusy Min/Max = %f / %f",((float)img.dimx())/((float)fmaxx)/2.0,((float)img.dimx())/((float)fminx)/2.0,((float)img.dimy())/((float)fmaxy)/2.0,((float)img.dimx())/((float)fminy)/2.0);
    hist.fill(0).draw_rectangle(0,0,thr,64,grey).draw_graph(visu[0].get_shared_channel(0).get_histogram(),white).draw_text(5,5,white,zero,11,0.6f,"Thr = %u",thr);
    visu[0].get_shared_channel(1) = 128*(1-visu[0].get_shared_channel(0).get_threshold(thr));
    visu[2].get_shared_channel(1) = (visu[0].get_shared_channel(0)-visu[0].get_shared_channel(0).get_dilate(3)).exp().threshold(0.9999).mul(visu[0].get_shared_channel(0).get_threshold(thr))*255;
    visu.display(disp);
    hist.display(histo);
  }
  return 0;
}
