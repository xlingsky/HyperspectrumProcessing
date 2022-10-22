#ifndef XLINGSKY_DESPECKLE_HPP
#define XLINGSKY_DESPECKLE_HPP

#include <queue>
#include <vector>

namespace xlingsky {
namespace raster {
namespace enhancement {
namespace detail {

template <typename T>
struct Histogram {
  typedef T value_type;
  value_type min;
  value_type step;
  int buckets;
  Histogram() : min(0), step(1), buckets(256) {}
  int bucket(value_type &v) const { 
      return (v-min)/step;
  }
};

typedef struct
{
  std::vector<int> elems;/* Number of pixels that fall into each luma bucket */
  std::queue<> origs;/* Original pixels */
  PixelsList origs[256]; 
  gint       xmin;
  gint       ymin;
  gint       xmax;
  gint       ymax; /* Source rect */
} DespeckleHistogram;

static inline void
list_add_elem (PixelsList   *list,
               const guchar *elem)
{
  const gint pos = list->start + list->count++;

  list->elems[pos >= MAX_LIST_ELEMS ? pos - MAX_LIST_ELEMS : pos] = elem;
}

static inline void
list_del_elem (PixelsList* list)
{
  list->count--;
  list->start++;

  if (list->start >= MAX_LIST_ELEMS)
    list->start = 0;
}

static inline const guchar *
list_get_random_elem (PixelsList *list)
{
  const gint pos = list->start + rand () % list->count;

  if (pos >= MAX_LIST_ELEMS)
    return list->elems[pos - MAX_LIST_ELEMS];

  return list->elems[pos];
}

static inline void
histogram_add (DespeckleHistogram *hist,
               guchar              val,
               const guchar       *orig)
{
  hist->elems[val]++;
  list_add_elem (&hist->origs[val], orig);
}

static inline void
histogram_remove (DespeckleHistogram *hist,
                  guchar              val)
{
  hist->elems[val]--;
  list_del_elem (&hist->origs[val]);
}

static inline void
histogram_clean (DespeckleHistogram *hist)
{
  gint i;

  for (i = 0; i < 256; i++)
    {
      hist->elems[i] = 0;
      hist->origs[i].count = 0;
    }
}

static inline const guchar *
histogram_get_median (DespeckleHistogram *hist,
                      const guchar       *_default)
{
  gint count = histrest;
  gint i;
  gint sum = 0;

  if (! count)
    return _default;

  count = (count + 1) / 2;

  i = 0;
  while ((sum += hist->elems[i]) < count)
    i++;

  return list_get_random_elem (&hist->origs[i]);
}

static inline void
add_val (DespeckleHistogram *hist,
         gint                black_level,
         gint                white_level,
         const guchar       *src,
         gint                width,
         gint                bpp,
         gint                x,
         gint                y)
{
  const gint pos   = (x + (y * width)) * bpp;
  const gint value = pixel_luminance (src + pos, bpp);

  if (value > black_level && value < white_level)
    {
      histogram_add (hist, value, src + pos);
      histrest++;
    }
  else
    {
      if (value <= black_level)
        hist0++;

      if (value >= white_level)
        hist255++;
    }
}

static inline void
del_val (DespeckleHistogram *hist,
         gint                black_level,
         gint                white_level,
         const guchar       *src,
         gint                width,
         gint                bpp,
         gint                x,
         gint                y)
{
  const gint pos   = (x + (y * width)) * bpp;
  const gint value = pixel_luminance (src + pos, bpp);

  if (value > black_level && value < white_level)
    {
      histogram_remove (hist, value);
      histrest--;
    }
  else
    {
      if (value <= black_level)
        hist0--;

      if (value >= white_level)
        hist255--;
    }
}

static inline void
add_vals (DespeckleHistogram *hist,
          gint                black_level,
          gint                white_level,
          const guchar       *src,
          gint                width,
          gint                bpp,
          gint                xmin,
          gint                ymin,
          gint                xmax,
          gint                ymax)
{
  gint x;
  gint y;

  if (xmin > xmax)
    return;

  for (y = ymin; y <= ymax; y++)
    {
      for (x = xmin; x <= xmax; x++)
        {
          add_val (hist,
                   black_level, white_level,
                   src, width, bpp, x, y);
        }
    }
}

static inline void
del_vals (DespeckleHistogram *hist,
          gint                black_level,
          gint                white_level,
          const guchar       *src,
          gint                width,
          gint                bpp,
          gint                xmin,
          gint                ymin,
          gint                xmax,
          gint                ymax)
{
  gint x;
  gint y;

  if (xmin > xmax)
    return;

  for (y = ymin; y <= ymax; y++)
    {
      for (x = xmin; x <= xmax; x++)
        {
          del_val (hist,
                   black_level, white_level,
                   src, width, bpp, x, y);
        }
    }
}

static inline void
update_histogram (DespeckleHistogram *hist,
                  gint                black_level,
                  gint                white_level,
                  const guchar       *src,
                  gint                width,
                  gint                bpp,
                  gint                xmin,
                  gint                ymin,
                  gint                xmax,
                  gint                ymax)
{
  /* assuming that radious of the box can change no more than one
     pixel in each call */
  /* assuming that box is moving either right or down */

  del_vals (hist,
            black_level, white_level,
            src, width, bpp, hist->xmin, hist->ymin, xmin - 1, hist->ymax);
  del_vals (hist,
            black_level, white_level,
            src, width, bpp, xmin, hist->ymin, xmax, ymin - 1);
  del_vals (hist,
            black_level, white_level,
            src, width, bpp, xmin, ymax + 1, xmax, hist->ymax);

  add_vals (hist,
            black_level, white_level,
            src, width, bpp, hist->xmax + 1, ymin, xmax, ymax);
  add_vals (hist,
            black_level, white_level,
            src, width, bpp, xmin, ymin, hist->xmax, hist->ymin - 1);
  add_vals (hist,
            black_level, white_level,
            src, width, bpp, hist->xmin, hist->ymax + 1, hist->xmax, ymax);

  hist->xmin = xmin;
  hist->ymin = ymin;
  hist->xmax = xmax;
  hist->ymax = ymax;
}
};
};  // namespace enhancement
};  // namespace raster
};  // namespace xlingsky

#endif