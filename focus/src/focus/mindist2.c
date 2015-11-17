#define xc100 TMx_xc[0]
#define yc100 TMx_xc[3]
#define zc100 TMx_xc[6]
#define xc010 TMx_xc[1]
#define yc010 TMx_xc[4]
#define zc010 TMx_xc[7]
#define xc001 TMx_xc[2]
#define yc001 TMx_xc[5]
#define zc001 TMx_xc[8]

  dx = dx0 - xc100; /* m00 */
  dy = dy0 - yc100;
  dz = dz0 - zc100;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc100; /* p00 */
  dy = dy0 + yc100;
  dz = dz0 + zc100;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc010; /* 0m0 */
  dy = dy0 - yc010;
  dz = dz0 - zc010;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc010; /* 0p0 */
  dy = dy0 + yc010;
  dz = dz0 + zc010;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc001; /* 00m */
  dy = dy0 - yc001;
  dz = dz0 - zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc001; /* 00p */
  dy = dy0 + yc001;
  dz = dz0 + zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc100 - xc010; /* mm0 */
  dy = dy0 - yc100 - yc010;
  dz = dz0 - zc100 - zc010;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc100 + xc010; /* pp0 */
  dy = dy0 + yc100 + yc010;
  dz = dz0 + zc100 + zc010;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc100 - xc001; /* m0m */
  dy = dy0 - yc100 - yc001;
  dz = dz0 - zc100 - zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc100 + xc001; /* p0p */
  dy = dy0 + yc100 + yc001;
  dz = dz0 + zc100 + zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc010 - xc001; /* 0mm */
  dy = dy0 - yc010 - yc001;
  dz = dz0 - zc010 - zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc010 + xc001; /* 0pp */
  dy = dy0 + yc010 + yc001;
  dz = dz0 + zc010 + zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc100 + xc010; /* mp0 */
  dy = dy0 - yc100 + yc010;
  dz = dz0 - zc100 + zc010;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc100 - xc010; /* pm0 */
  dy = dy0 + yc100 - yc010;
  dz = dz0 + zc100 - zc010;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc100 + xc001; /* m0p */
  dy = dy0 - yc100 + yc001;
  dz = dz0 - zc100 + zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc100 - xc001; /* p0m */
  dy = dy0 + yc100 - yc001;
  dz = dz0 + zc100 - zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc010 + xc001; /* 0mp */
  dy = dy0 - yc010 + yc001;
  dz = dz0 - zc010 + zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc010 - xc001; /* 0pm */
  dy = dy0 + yc010 - yc001;
  dz = dz0 + zc010 - zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc100 - xc010 - xc001; /* mmm */
  dy = dy0 - yc100 - yc010 - yc001;
  dz = dz0 - zc100 - zc010 - zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc100 + xc010 + xc001; /* ppp */
  dy = dy0 + yc100 + yc010 + yc001;
  dz = dz0 + zc100 + zc010 + zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc100 - xc010 + xc001; /* mmp */
  dy = dy0 - yc100 - yc010 + yc001;
  dz = dz0 - zc100 - zc010 + zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc100 + xc010 - xc001; /* ppm */
  dy = dy0 + yc100 + yc010 - yc001;
  dz = dz0 + zc100 + zc010 - zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc100 + xc010 - xc001; /* mpm */
  dy = dy0 - yc100 + yc010 - yc001;
  dz = dz0 - zc100 + zc010 - zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc100 - xc010 + xc001; /* pmp */
  dy = dy0 + yc100 - yc010 + yc001;
  dz = dz0 + zc100 - zc010 + zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 - xc100 + xc010 + xc001; /* mpp */
  dy = dy0 - yc100 + yc010 + yc001;
  dz = dz0 - zc100 + zc010 + zc001;
  EvalDist2(dx, dy, dz);

  dx = dx0 + xc100 - xc010 - xc001; /* pmm */
  dy = dy0 + yc100 - yc010 - yc001;
  dz = dz0 + zc100 - zc010 - zc001;
  EvalDist2(dx, dy, dz);

#undef xc100
#undef yc100
#undef zc100
#undef xc010
#undef yc010
#undef zc010
#undef xc001
#undef yc001
#undef zc001
