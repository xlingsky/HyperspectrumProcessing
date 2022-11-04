
namespace xlingsky{
namespace raster{

template<class Operator>
#ifdef TRANSFORM_OMP
void transform_omp
#else
void transform
#endif
(
#ifdef TRANSFORM_CONST
    const void* data,
#else
    void* data,
#endif
    int cols, int rows, int colspace, int rowspace, Operator& op){
#ifdef TRANSFORM_CONST
  const char* pdata = (const char*)data;
#else
  char* pdata = (char*)data;
#endif

#ifdef TRANSFORM_OMP
#pragma omp parallel for
#endif
  for(int r=0; r<rows; ++r){
#ifdef TRANSFORM_CONST
    const char* pr = pdata+r*rowspace;
#else
    char* pr = pdata+r*rowspace;
#endif
    for(int c=0; c<cols; ++c){
      op(pr+c*colspace);
    }
  }
}

};
};
