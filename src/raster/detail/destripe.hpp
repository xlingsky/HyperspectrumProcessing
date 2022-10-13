#ifndef XLINGSKY_DESTRIPE_HPP
#define XLINGSKY_DESTRIPE_HPP

bool Destripe(int tile_width){
  GeglBuffer *src_buffer;
  GeglBuffer *dest_buffer;
  const Babl *format;
  guchar     *src_rows;       /* image data */
  gint        x1, x2, y1;
  gint        width, height;
  gint        bpp;
  glong      *hist, *corr;        /* "histogram" data */
  gint        tile_width = gimp_tile_width ();
  gint        i, x, y, ox, cols;

  x2 = x1 + width;

  if (gimp_drawable_is_rgb (drawable))
    {
      if (gimp_drawable_has_alpha (drawable))
        format = babl_format ("R'G'B'A u8");
      else
        format = babl_format ("R'G'B' u8");
    }
  else
    {
      if (gimp_drawable_has_alpha (drawable))
        format = babl_format ("Y'A u8");
      else
        format = babl_format ("Y' u8");
    }

  bpp = babl_format_get_bytes_per_pixel (format);

  /*
   * Setup for filter...
   */

  src_buffer  = gimp_drawable_get_buffer (drawable);
  dest_buffer = gimp_drawable_get_shadow_buffer (drawable);

  hist = g_new (long, width * bpp);
  corr = g_new (long, width * bpp);
  src_rows = g_new (guchar, tile_width * height * bpp);

  memset (hist, 0, width * bpp * sizeof (long));

  /*
   * collect "histogram" data.
   */

  for (ox = x1; ox < x2; ox += tile_width)
    {
      guchar *rows = src_rows;

      cols = x2 - ox;
      if (cols > tile_width)
        cols = tile_width;

      gegl_buffer_get (src_buffer, GEGL_RECTANGLE (ox, y1, cols, height), 1.0,
                       format, rows,
                       GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);

      for (y = 0; y < height; y++)
        {
          long   *h       = hist + (ox - x1) * bpp;
          guchar *row_end = rows + cols * bpp;

          while (rows < row_end)
            *h++ += *rows++;
        }

      if (! preview)
        gimp_progress_update (progress += progress_inc);
    }

  /*
   * average out histogram
   */

  {
    gint extend = (vals.avg_width / 2) * bpp;

    for (i = 0; i < MIN (3, bpp); i++)
      {
        long *h   = hist - extend + i;
        long *c   = corr - extend + i;
        long  sum = 0;
        gint  cnt = 0;

        for (x = -extend; x < width * bpp; x += bpp)
          {
            if (x + extend < width * bpp)
              {
                sum += h[ extend]; cnt++;
              }

            if (x - extend >= 0)
              {
                sum -= h[-extend]; cnt--;
              }

            if (x >= 0)
              {
                if (*h)
                  *c = ((sum / cnt - *h) << 10) / *h;
                else
                  *c = G_MAXINT;
              }

            h += bpp;
            c += bpp;
          }
      }
  }

  /*
   * remove stripes.
   */

  for (ox = x1; ox < x2; ox += tile_width)
    {
      guchar *rows = src_rows;

      cols = x2 - ox;
      if (cols > tile_width)
        cols = tile_width;

      gegl_buffer_get (src_buffer, GEGL_RECTANGLE (ox, y1, cols, height), 1.0,
                       format, rows,
                       GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);

      if (! preview)
        gimp_progress_update (progress += progress_inc);

      for (y = 0; y < height; y++)
        {
          long   *c = corr + (ox - x1) * bpp;
          guchar *row_end = rows + cols * bpp;

          if (vals.histogram)
            {
              while (rows < row_end)
                {
                  *rows = MIN (255, MAX (0, 128 + (*rows * *c >> 10)));
                  c++; rows++;
                }
            }
          else
            {
              while (rows < row_end)
                {
                  *rows = MIN (255, MAX (0, *rows + (*rows * *c >> 10) ));
                  c++; rows++;
                }
            }
        }

      gegl_buffer_set (dest_buffer, GEGL_RECTANGLE (ox, y1, cols, height), 0,
                       format, src_rows,
                       GEGL_AUTO_ROWSTRIDE);

      if (! preview)
        gimp_progress_update (progress += progress_inc);
    }

  g_free (src_rows);

  g_object_unref (src_buffer);

  if (preview)
    {
      guchar *buffer = g_new (guchar, width * height * bpp);

      gegl_buffer_get (dest_buffer, GEGL_RECTANGLE (x1, y1, width, height), 1.0,
                       format, buffer,
                       GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);

      gimp_preview_draw_buffer (GIMP_PREVIEW (preview),
                                buffer, width * bpp);

      g_free (buffer);
      g_object_unref (dest_buffer);
    }
  else
    {
      g_object_unref (dest_buffer);

      gimp_progress_update (1.0);

      gimp_drawable_merge_shadow (drawable, TRUE);
      gimp_drawable_update (drawable,
                            x1, y1, width, height);
    }

  g_free (hist);
  g_free (corr);
}

#endif

