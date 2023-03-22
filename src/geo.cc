#include "pos.h"

#include "pathadaptor.hpp"

int main(int argc, char *argv[]) {
  boost::filesystem::path path(argv[1]);
  HSP::Pos pos;
  pos.load(argv[1]);
  HSP::PinholeCamera cam;
  cam.SetFocalLength(59750/30.0);
  cam.SetPrincipalPoint(-724,0);
  double a[3] = {0,0,0};
  cam.SetAngles(a);
  double t[3] = {0.006, 0.25768, -0.6773 };
  cam.SetTranslation(t);
  HSP::LinescanModel model;
  model.SetCamera(&cam);
  model.SetPos(&pos);
  boost::filesystem::path rpcpath(path);
  rpcpath.replace_extension(".rpc");
  model.test();
  double range_samp[] = {10, 1300};
  double range_line[] = {8, (double)pos.data().rbegin()->first-8};
  double range_height[] = {1500, 2500};
  model.GenerateRPC(range_samp, range_line, range_height, path.string().c_str());
  return 0;
}
