#ifndef XLINGSKY_DESTRIPE_HPP
#define XLINGSKY_DESTRIPE_HPP

#include <vector>
#include <limits>
#include <type_traits>

/*
* destripe along column
*/
template<typename DataType, typename HistType = long>
bool Destripe(
    int cols, int rows,
    void* src, int src_col_space, int src_row_space,
    void* dst, int dst_col_space, int dst_row_space, int win_width,
    HistType hist_mean = std::numeric_limits<typename HistType>::min()
){
  std::vector<HistType> hist(cols, 0), corr(cols);

  char *psrc, *pdst;
  /*
   * collect "histogram" data.
   */
  psrc = (char*)src;
  for (int c=0; c < cols; ++c) {
    char *pr_src = psrc;
    for (int r = 0; r < rows; r++) {
      hist[c] += (HistType)*((DataType*)pr_src);
      pr_src += src_row_space;
    }
    psrc += src_col_space;
   }

  /*
   * average out histogram
   */

  {
     int extend = win_width / 2;
     HistType *h = &hist[0] - extend;
     HistType *c = &corr[0] - extend;
     HistType sum = 0;
     int cnt = 0;

     for (int x = -extend; x < cols; ++x) {
       if (x + extend < cols) {
         sum += h[extend];
         cnt++;
       }

       if (x - extend >= 0) {
         sum -= h[-extend];
         cnt--;
       }

       if (x >= 0) {
         if (*h) {
           if (
               std::is_same<int, HistType>::value || std::is_same<unsigned int, HistType>::value ||
               std::is_same<long, HistType>::value || std::is_same<unsigned long, HistType>::value ||
               std::is_same<size_t, HistType>::value
               )
             *c = ((sum / cnt - *h) << 10) / *h;
           else
             *c = (sum / cnt - *h) / *h;
         }
         else
           *c = std::numeric_limits<HistType>::max();
       }

       ++h;
       ++c;
     }
  }

  /*
   * remove stripes.
   */
  psrc = (char*)src;
  pdst = (char*)dst;
  for (int x=0; x < cols; ++x) {
    HistType &c = corr[x];
    char *pr_src = psrc;
    char *pr_dst = pdst;
    if (hist_mean==std::numeric_limits<typename HistType>::min()) {
      for (int y = 0; y < rows; ++y) {
        DataType &v = *((DataType *)pr_src);
        if (std::is_same<int, HistType>::value ||
            std::is_same<unsigned int, HistType>::value ||
            std::is_same<long, HistType>::value ||
            std::is_same<unsigned long, HistType>::value ||
            std::is_same<size_t, HistType>::value) {
          *((DataType *)pr_dst) =
              (DataType)(std::min)(std::numeric_limits<DataType>::max(), std::max(0,v+(v*c>>10)) );
        } else {
          *((DataType *)pr_dst) =
              (DataType)(std::min)(std::numeric_limits<DataType>::max(), std::max(0,v+(v*c)) );
        }
        pr_src += src_row_space;
        pr_dst += dst_row_space;
      }
    } else {
      for (int y = 0; y < rows; ++y) {
        DataType &v = *((DataType *)pr_src);
        if (std::is_same<int, HistType>::value ||
            std::is_same<unsigned int, HistType>::value ||
            std::is_same<long, HistType>::value ||
            std::is_same<unsigned long, HistType>::value ||
            std::is_same<size_t, HistType>::value) {
          *((DataType *)pr_dst) =
              (DataType)(std::min)(std::numeric_limits<DataType>::max(), std::max(0,hist_mean+(v*c>>10)) );
        } else {
          *((DataType *)pr_dst) =
              (DataType)(std::min)(std::numeric_limits<DataType>::max(), std::max(0,hist_mean+(v*c)) );
        }
        pr_src += src_row_space;
        pr_dst += dst_row_space;
      }
    }
    psrc += src_col_space;
    pdst += dst_col_space;
  }

  return true;
}

#endif

