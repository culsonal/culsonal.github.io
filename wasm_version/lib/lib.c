
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include <stdio.h>

typedef struct { float x; float y; } vec2;
typedef struct { float x; float y; float z; } vec3;
typedef struct { float m11; float m12; float m21; float m22; } mat2x2;
typedef struct {
  float m11; float m12; float m13;
  float m21; float m22; float m23;
  float m31; float m32; float m33;
} mat3x3;

float*    alloc_float32array(int size) { return (float*)calloc(size, sizeof(float)); }
uint8_t*  alloc_uint8array(int size)   { return (uint8_t*)calloc(size, sizeof(uint8_t)); }
uint16_t* alloc_uint16array(int size)  { return (uint16_t*)calloc(size, sizeof(uint16_t)); }
uint32_t* alloc_uint32array(int size)  { return (uint32_t*)calloc(size, sizeof(uint32_t)); }

inline float rnd_uniform() { return (rand() % 10000)/10000.; }
inline float rnd() { return rnd_uniform()*2.-1.; }

inline float clampf(float v, float min, float max) { return fmax(min, fmin(v, max)); }
inline float dot2(vec2 v, vec2 w) { return v.x*w.x + v.y*w.y; }
inline float dot3(vec3 v, vec3 w) { return v.x*w.x + v.y*w.y + v.z*w.z; }
inline vec2 add2(vec2 v, vec2 w)  { return (vec2){ v.x+w.x, v.y+w.y }; }
inline vec3 add3(vec3 v, vec3 w)  { return (vec3){ v.x+w.x, v.y+w.y, v.z+w.z }; }
inline vec2 sub2(vec2 v, vec2 w)  { return (vec2){ v.x-w.x, v.y-w.y }; }
inline vec3 sub3(vec3 v, vec3 w)  { return (vec3){ v.x-w.x, v.y-w.y, v.z-w.z }; }
inline vec3 cross(vec3 v, vec3 w) { return (vec3){ v.y*w.z-v.z*w.y, v.z*w.x-v.x*w.z, v.x*w.y-v.y*w.x }; }
inline float len2(vec2 v) { return sqrtf(v.x*v.x + v.y*v.y); }
inline float len3(vec3 v) { return sqrtf(v.x*v.x + v.y*v.y + v.z*v.z); }
inline float dist2(vec2 v, vec2 w) { return len2(sub2(v, w)); }
inline float dist3(vec3 v, vec3 w) { return len3(sub3(v, w)); }

inline float det2(mat2x2 m) { return m.m11*m.m22 - m.m12*m.m21; }
inline mat2x2 mat_mat_prod2x2(mat2x2 m1, mat2x2 m2) { return (mat2x2){
  m1.m11*m2.m11 + m1.m12*m2.m21, m1.m11*m2.m12 + m1.m12*m2.m22,
  m1.m21*m2.m11 + m1.m22*m2.m21, m1.m21*m2.m12 + m1.m22*m2.m22 }; }
inline float tri_area_2d_signed(float x1, float y1, float x2, float y2, float x3, float y3)
  { return (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1); }

inline vec3 rot_xz(vec3 v, float ang)
  { return (vec3){ cosf(ang)*v.x - sinf(ang)*v.z,   v.y,   sinf(ang)*v.x + cosf(ang)*v.z }; }
inline vec3 rot_yz(vec3 v, float ang)
  { return (vec3){ v.x,   cosf(ang)*v.y - sinf(ang)*v.z,   sinf(ang)*v.y + cosf(ang)*v.z }; }
inline vec3 rot_xy(vec3 v, float ang)
  { return (vec3){ cosf(ang)*v.x - sinf(ang)*v.y,   sinf(ang)*v.x + cosf(ang)*v.y,   v.z }; }

void fill_arr_pt_cloud_3d_sphere(float* pts, int num_pts, float radius) { /*todo*/ }
void fill_arr_pt_cloud_2d_cube(float* pts, int num_pts, float cube_size)
  { srand(3); for (int i = 0; i < num_pts; i++) for (int j = 0; j < 2; j++) pts[i*2+j] = rnd()*cube_size; }
void fill_arr_pt_cloud_3d_cube(float* pts, int num_pts, float cube_size)
  { srand(3); for (int i = 0; i < num_pts; i++) for (int j = 0; j < 3; j++) pts[i*3+j] = rnd()*cube_size; }

void move_pts_2d(float* pts, int num_pts, float ox, float oy)
  { for (int i = 0; i < num_pts; i++) { pts[i*2+0] += ox; pts[i*2+1] += oy; } }
void move_pts_3d(float* pts, int num_pts, float ox, float oy, float oz)
  { for (int i = 0; i < num_pts; i++) { pts[i*3+0] += ox; pts[i*3+1] += oy; pts[i*3+2] += oz; } }

// I think it implements curl noise right but it doesn't look quite windy, but I'll let it stay
void modify_vels_via_curl_noise_2d(float* pts, int num_pts, float* inv_mass, float* vels, float dt, int num_basis_fns,
    float strength, int noise_seed) {
  srand(noise_seed);
  for (int i = 0; i < num_pts; i++) {
    const float x = pts[i*2+0];
    const float y = pts[i*2+1];
    float dx = 0; float dy = 0;
    for (int j = 0; j < num_basis_fns; j++) {
      const float freq1   = 20.*rnd();     const float freq2   = 20.*rnd();
      const float phase1  = 3.1415*rnd(); const float phase2  = 3.1415*rnd();
      const float weight = rnd()*.1;
      dx += freq1*weight*cosf(freq1*x + phase1)*sinf(freq2*y + phase2);
      dy += freq2*weight*cosf(freq2*y + phase2)*sinf(freq1*x + phase1);
    }

    const float l = sqrtf(dx*dx + dy*dy);
    const float orth_dx = -dy / l;
    const float orth_dy =  dx / l;
    vels[i*2+0] += orth_dx*strength*dt*inv_mass[i];
    vels[i*2+1] += orth_dy*strength*dt*inv_mass[i];
  }
}

inline int hash_fn_pt_2d(float x, float y, int hash_table_size, int dim, float bucket_size) {
  const int xi = (int)floor((x*.5+.5)/bucket_size);
  const int yi = (int)floor((y*.5+.5)/bucket_size);
  return llabs((xi * 92837111) ^ (yi * 689287499)) % hash_table_size;
}

inline int hash_fn_pt_3d(float x, float y, float z, int hash_table_size, int dim, float bucket_size) {
  const int xi = (int)floor((x*.5+.5)/bucket_size);
  const int yi = (int)floor((y*.5+.5)/bucket_size);
  const int zi = (int)floor((z*.5+.5)/bucket_size);
  return llabs((xi * 92837111) ^ (yi * 689287499) ^ (zi * 283923481)) % hash_table_size;
}

inline int hash_fn(float* pts, int idx, int hash_table_size, int dim, float bucket_size) {
  return dim == 2 ?
    hash_fn_pt_2d(pts[idx*dim+0], pts[idx*dim+1], hash_table_size, dim, bucket_size) :
    hash_fn_pt_3d(pts[idx*dim+0], pts[idx*dim+1], pts[idx*dim+2], hash_table_size, dim, bucket_size);
}

void compute_hash(float* pts, int num_pts, uint16_t* pts_reordered, uint16_t* hash, float particle_radius, int dim) {
  const int num_elems_in_bucket = 0;
  const int bucket_filled_to_idx = 1;
  const int sum_of_previous_bucket_counts = 2;
  const int hash_table_size = num_pts; // assumption

  for (int i = 0; i < hash_table_size*3; i++) hash[i] = 0;
  for (int i = 0; i < num_pts; i++) {
    int h = hash_fn(pts, i, hash_table_size, dim, particle_radius);
    hash[h*3+num_elems_in_bucket] += 1;
  }

  int sum = 0;
  for (int i = 0; i < hash_table_size; i++) {
    hash[i*3+sum_of_previous_bucket_counts] = sum;
    sum += hash[i*3+num_elems_in_bucket];
  }

  for (int i = 0; i < num_pts; i++) {
    int h = hash_fn(pts, i, hash_table_size, dim, particle_radius);
    int idx_reordered = hash[h*3+sum_of_previous_bucket_counts] + hash[h*3+bucket_filled_to_idx];
    pts_reordered[idx_reordered] = i;
    hash[h*3+bucket_filled_to_idx] += 1;
  }
}

// todo offset, aspect ratio
void draw_pts_2d_with_colors(float* pts, int num_pts, uint8_t* img, int img_width, int img_height,
    int pt_size, float* colors) {
  const int dim = 2;
  for (int i = 0; i < img_width*img_height*4; i++) img[i] = 0;
  for (int i = 0; i < num_pts; i++) {
    const int x = (int)((pts[i*dim+0]*.5+.5)*(float)img_width);
    const int y = (int)((pts[i*dim+1]*.5+.5)*(float)img_height);
    for (int dx = -(pt_size-1); dx <= (pt_size-1); dx++) for (int dy = -(pt_size-1); dy <= (pt_size-1); dy++) {
      const int pix_idx = (y+dy)*img_height + (x+dx);
      if (pix_idx < 0 || pix_idx > img_width*img_height*4) continue;
      img[pix_idx*4+0] = (uint8_t)(colors[i*3+0]*255);
      img[pix_idx*4+1] = (uint8_t)(colors[i*3+1]*255);
      img[pix_idx*4+2] = (uint8_t)(colors[i*3+2]*255);
      img[pix_idx*4+3] = 255;
    }
  }
}

void draw_pts_2d(float* pts, int num_pts, uint8_t* img, int img_width, int img_height, int pt_size) {
  float colors[num_pts*3];
  srand(5);
  for (int i = 0; i < num_pts*3; i++) colors[i] = rnd_uniform();
  draw_pts_2d_with_colors(pts, num_pts, img, img_width, img_height, pt_size, colors);
}

void draw_surface_pts_2d(float* pts, int num_pts, uint16_t* pts_reordered, uint16_t* hash, float particle_radius,
  uint8_t* img_ptr, int img_width, int img_height, int pt_size) { /* gutted; use the proper illum algo maybe */ }

void draw_pts_3d(float* pts, int num_pts, uint8_t* img, float* depth, int img_width, int img_height, int pt_size,
    float offx, float offy, float offz, float rot_xz_ang, float rot_xy_ang, float rot_yz_ang, float max_depth) {
  const int dim = 3;
  for (int i = 0; i < img_width*img_height*4; i++) img[i] = 0;
  for (int i = 0; i < img_width*img_height; i++) depth[i] = 0;

  // depth
  for (int i = 0; i < num_pts; i++) {
    //vec3 pt = (vec3){ pts[i*dim+0]+offx, pts[i*dim+1]+offy, pts[i*dim+2]+offz };
    vec3 pt = (vec3){ pts[i*dim+0], pts[i*dim+1], pts[i*dim+2] };
    pt = rot_xz(pt, rot_xz_ang);
    pt = rot_xy(pt, rot_xy_ang);
    pt = rot_yz(pt, rot_yz_ang);
    pt = (vec3){ pt.x+offx, pt.y+offy, pt.z+offz };
    const int x = (int)((pt.x/pt.z*.5+.5)*(float)img_width);
    const int y = (int)((pt.y/pt.z*.5+.5)*(float)img_height);
    const float depth_val = 1.-pt.z/max_depth;
    if (pt.z < 0) continue;
    int pt_size_aug = (int)(pt_size*powf(depth_val, 4));
    for (int dx = -pt_size_aug; dx <= pt_size_aug; dx++) for (int dy = -pt_size_aug; dy <= pt_size_aug; dy++) {
      const int pix_idx = (y+dy)*img_width + (x+dx);
      if (pix_idx < 0 || pix_idx > img_width*img_height) continue;
      depth[pix_idx] = fmax(depth[pix_idx], depth_val);
    }
  }

  // rgb
  srand(5);
  for (int i = 0; i < num_pts; i++) {
    //vec3 pt = (vec3){ pts[i*dim+0]+offx, pts[i*dim+1]+offy, pts[i*dim+2]+offz };
    vec3 pt = (vec3){ pts[i*dim+0], pts[i*dim+1], pts[i*dim+2] };
    pt = rot_xz(pt, rot_xz_ang);
    pt = rot_xy(pt, rot_xy_ang);
    pt = rot_yz(pt, rot_yz_ang);
    pt = (vec3){ pt.x+offx, pt.y+offy, pt.z+offz };
    const int x = (int)((pt.x/pt.z*.5+.5)*(float)img_width);
    const int y = (int)((pt.y/pt.z*.5+.5)*(float)img_height);
    const float depth_val = 1.-pt.z/max_depth;
    vec3 col = (vec3){ .6+.4*rnd_uniform(), .6+.4*rnd_uniform(), .6+.4*rnd_uniform() };
    int pt_size_aug = (int)(pt_size*powf(depth_val, 4));
    for (int dx = -pt_size_aug; dx <= pt_size_aug; dx++) for (int dy = -pt_size_aug; dy <= pt_size_aug; dy++) {
      const int pix_idx = (y+dy)*img_width + (x+dx);
      if (pix_idx < 0 || pix_idx >= img_width*img_height*4) continue;
      if (depth_val < 0 || depth[pix_idx] > depth_val) continue;

      float shade = powf(depth_val, 2.);
      vec3 shaded_col = (vec3){ col.x*shade, col.y*shade, col.z*shade };
      img[pix_idx*4+0] = (uint8_t)(shaded_col.x * 255.);
      img[pix_idx*4+1] = (uint8_t)(shaded_col.y * 255.);
      img[pix_idx*4+2] = (uint8_t)(shaded_col.z * 255.);
      img[pix_idx*4+3] = 255;
    }
  }
}

void draw_pts_3d_ray_traced(float* pts, int num_pts, uint8_t* img, int img_width, int img_height, int pt_size,
  float offx, float offy, float offz, float rot_xz_ang, float rot_xy_ang) { /* todo */ }

void pbd_pre(float* pts, int num_pts, float* prev, float* vels, float dt, int num_substeps, int dim) {
  const float sub_frame_dt = dt/(float)(num_substeps);
  for (int i = 0; i < num_pts*dim; i++) prev[i] = pts[i];
  for (int i = 0; i < num_pts*dim; i++) pts[i] += sub_frame_dt*vels[i];
}

void pbd_post(float* pts, int num_pts, float* prev, float* vels, float* orig, float* inv_mass, float dt,
    int num_substeps, int dim, float grav, float particle_radius, float air_mass, float vel_damping, float clamp_val,
    bool clamp_box, bool limit_max_vel) {
  const float sub_frame_dt = dt/(float)(num_substeps);
  if (clamp_box) for (int i = 0; i < num_pts; i++) pts[i*dim+0] = clampf(pts[i*dim+0], -clamp_val, clamp_val);
  if (clamp_box) for (int i = 0; i < num_pts; i++) pts[i*dim+1] = clampf(pts[i*dim+1], -clamp_val, clamp_val);
  if (clamp_box) for (int i = 0; i < num_pts; i++) pts[i*dim+2] = clampf(pts[i*dim+2], -clamp_val, clamp_val);

  for (int i = 0; i < num_pts*dim; i++) vels[i] = (pts[i] - prev[i])/sub_frame_dt;
  for (int i = 0; i < num_pts; i++) vels[i*dim+1] += (inv_mass[i]-air_mass)*grav;
  for (int i = 0; i < num_pts*dim; i++) vels[i] *= vel_damping;
  if (limit_max_vel) for (int i = 0; i < num_pts; i++) {
    const float max_vel = .5*particle_radius/sub_frame_dt;
    float vel_strength = 0.;
    if (dim == 2) { vec2 v = (vec2){ vels[i*dim+0], vels[i*dim+1] }; vel_strength = len2(v); }
    if (dim == 3) { vec3 v = (vec3){ vels[i*dim+0], vels[i*dim+1], vels[i*dim+2] }; vel_strength = len3(v); }
    if (vel_strength > max_vel) for (int j = 0; j < 3; j++) vels[i*dim+j] *= max_vel/vel_strength;
  }
}

void solve_length_constraints_of_pt_pairs_2d(float* pts, float* inv_mass, int i, int j, float desired_dist,
    float damping) {
  const float x1 = pts[i*2+0]; const float y1 = pts[i*2+1];
  const float x2 = pts[j*2+0]; const float y2 = pts[j*2+1];
  const float im1 = inv_mass[i]; const float im2 = inv_mass[j];
  const vec2 p1 = (vec2){ x1, y1 }; const vec2 p2 = (vec2){ x2, y2 };
  const float dist = dist2(p1, p2);
  const float d1x = (x1-x2)/dist; const float d2x = -d1x;
  const float d1y = (y1-y2)/dist; const float d2y = -d1y;
  const float grad_sum_sq = im1*(d1x*d1x + d1y*d1y) + im2*(d2x*d2x + d2y*d2y);
  if (grad_sum_sq < .0000001) return;
  const float constraint = dist - desired_dist;
  const float lambda = -constraint*damping/grad_sum_sq;
  pts[i*2+0] += im1*lambda*d1x; pts[i*2+1] += im1*lambda*d1y;
  pts[j*2+0] += im2*lambda*d2x; pts[j*2+1] += im2*lambda*d2y;
}

void solve_length_constraints_of_pt_pairs_3d(float* pts, float* inv_mass, int i, int j, float desired_dist,
    float damping) {
  const float x1 = pts[i*3+0]; const float y1 = pts[i*3+1]; const float z1 = pts[i*3+2];
  const float x2 = pts[j*3+0]; const float y2 = pts[j*3+1]; const float z2 = pts[j*3+2];
  const float im1 = inv_mass[i]; const float im2 = inv_mass[j];
  const vec3 p1 = (vec3){ x1, y1, z1 }; const vec3 p2 = (vec3){ x2, y2, z2 };
  const float dist = dist3(p1, p2);
  const float d1x = (x1-x2)/dist; const float d2x = -d1x;
  const float d1y = (y1-y2)/dist; const float d2y = -d1y;
  const float d1z = (z1-z2)/dist; const float d2z = -d1z;
  const float grad_sum_sq = im1*(d1x*d1x + d1y*d1y + d1z*d1z) + im2*(d2x*d2x + d2y*d2y + d2z*d2z);
  if (grad_sum_sq < .0000001) return;
  const float constraint = dist - desired_dist;
  const float lambda = -constraint*damping/grad_sum_sq;
  pts[i*3+0] += im1*lambda*d1x; pts[i*3+1] += im1*lambda*d1y; pts[i*3+2] += im1*lambda*d1z;
  pts[j*3+0] += im2*lambda*d2x; pts[j*3+1] += im2*lambda*d2y; pts[j*3+2] += im2*lambda*d2z;
}

void solve_pairwise_constraints_2d(float* pts, int num_pts, uint16_t* pts_reordered, float* inv_mass, uint16_t* hash,
    float particle_radius, float damping, float sticky_mul) {
  const int num_elems_in_bucket = 0;
  const int bucket_filled_to_idx = 1;
  const int sum_of_previous_bucket_counts = 2;
  const int dim = 2;
  for (int i = 0; i < num_pts; i++) {
    const float x = pts[i*dim+0];
    const float y = pts[i*dim+1];
    for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) {
      const int h = hash_fn_pt_2d(x+dx*particle_radius*2, y+dy*particle_radius*2, num_pts, dim, particle_radius);
      const int bucket_start_idx = hash[h*3+sum_of_previous_bucket_counts];
      const int bucket_end_idx = bucket_start_idx + hash[h*3+num_elems_in_bucket];
      for (int j = bucket_start_idx; j < bucket_end_idx; j++) {
        const int neigh_idx = pts_reordered[j];
        if (i == neigh_idx) continue;
        const vec2 p = (vec2){ x, y }; const vec2 np = (vec2){ pts[neigh_idx*dim+0], pts[neigh_idx*dim+1] };
        if (dist2(p, np) < particle_radius)
          solve_length_constraints_of_pt_pairs_2d(pts, inv_mass, i, neigh_idx, particle_radius*sticky_mul, damping);
  }}}
}

void solve_pairwise_constraints_3d(float* pts, int num_pts, uint16_t* pts_reordered, float* inv_mass, uint16_t* hash,
    float particle_radius, float damping, float sticky_mul) {
  const int num_elems_in_bucket = 0;
  const int bucket_filled_to_idx = 1;
  const int sum_of_previous_bucket_counts = 2;
  const int dim = 3;
  for (int i = 0; i < num_pts; i++) {
    const float x = pts[i*dim+0];
    const float y = pts[i*dim+1];
    const float z = pts[i*dim+2];
    for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) for (int dz = -1; dz <= 1; dz++) {
      const int h = hash_fn_pt_3d(
        x+dx*particle_radius*2,
        y+dy*particle_radius*2,
        z+dz*particle_radius*2, num_pts, dim, particle_radius);
      const int bucket_start_idx = hash[h*3+sum_of_previous_bucket_counts];
      const int bucket_end_idx = bucket_start_idx + hash[h*3+num_elems_in_bucket];
      for (int j = bucket_start_idx; j < bucket_end_idx; j++) {
        const int neigh_idx = pts_reordered[j];
        if (i == neigh_idx) continue;
        const vec3  p = (vec3){ x, y, z };
        const vec3 np = (vec3){ pts[neigh_idx*dim+0], pts[neigh_idx*dim+1], pts[neigh_idx*dim+2] };
        if (dist3(p, np) < particle_radius)
          solve_length_constraints_of_pt_pairs_3d(pts, inv_mass, i, neigh_idx, particle_radius*sticky_mul, damping);
  }}}
}

void solve_neohookean_constraints_tris_2d(float* pts, int num_pts, uint16_t* tris, int num_tris, float* orig,
    uint8_t* is_tri_active, float* inv_mass, float ch_damping, float cd_damping) {
  const int dim = 2;
  for (int i = 0; i < num_tris; i++) {
    if (is_tri_active[i] == 0) continue;
    const int i1 = tris[i*3+0]; const int i2 = tris[i*3+1]; const int i3 = tris[i*3+2];
    const float  x1 = pts[i1*dim+0];  const float  x2 = pts[i2*dim+0];  const float  x3 = pts[i3*dim+0];
    const float  y1 = pts[i1*dim+1];  const float  y2 = pts[i2*dim+1];  const float  y3 = pts[i3*dim+1];
    const float ox1 = orig[i1*dim+0]; const float ox2 = orig[i2*dim+0]; const float ox3 = orig[i3*dim+0];
    const float oy1 = orig[i1*dim+1]; const float oy2 = orig[i2*dim+1]; const float oy3 = orig[i3*dim+1];
    const float im1 = inv_mass[i1];   const float im2 = inv_mass[i2];   const float im3 = inv_mass[i3];
    mat2x2 X    = (mat2x2){   x2-x1,   x3-x1,   y2-y1,   y3-y1 };
    mat2x2 Xhat = (mat2x2){ ox2-ox1, ox3-ox1, oy2-oy1, oy3-oy1 };
    const float detB = 1./det2(Xhat);
    const float b1 =  detB*Xhat.m22; const float b2 = -detB*Xhat.m12;
    const float b3 = -detB*Xhat.m21; const float b4 =  detB*Xhat.m11;
    const mat2x2 B = (mat2x2){ b1, b2, b3, b4 };
    const mat2x2 F = mat_mat_prod2x2(X, B);
    const float dCHdX1x = detB*(y2 - y3); const float dCHdX1y = detB*(x3 - x2);
    const float dCHdX2x = detB*(y3 - y1); const float dCHdX2y = detB*(x1 - x3);
    const float dCHdX3x = detB*(y1 - y2); const float dCHdX3y = detB*(x2 - x1);
    const float ch_grad_sum_sq =
      im1*(dCHdX1x*dCHdX1x + dCHdX1y*dCHdX1y) +
      im2*(dCHdX2x*dCHdX2x + dCHdX2y*dCHdX2y) +
      im3*(dCHdX3x*dCHdX3x + dCHdX3y*dCHdX3y);
    const float ch_constraint = det2(F) - (1.+cd_damping/ch_damping);
    //const float ch_constraint = det2(F) - 1.;
    const float ch_pbd_lambda = -ch_constraint*ch_damping/ch_grad_sum_sq;
    const float dx1h = im1*ch_pbd_lambda*dCHdX1x; const float dy1h = im1*ch_pbd_lambda*dCHdX1y;
    const float dx2h = im2*ch_pbd_lambda*dCHdX2x; const float dy2h = im2*ch_pbd_lambda*dCHdX2y;
    const float dx3h = im3*ch_pbd_lambda*dCHdX3x; const float dy3h = im3*ch_pbd_lambda*dCHdX3y;

    const float f1 = F.m11; const float f2 = F.m12;
    const float f3 = F.m21; const float f4 = F.m22;
    const float cd_constraint = sqrtf(f1*f1 + f2*f2 + f3*f3 + f4*f4);
    const float dCDdX1x = -(f1*b1 + f1*b3 + f2*b2 + f2*b4)/cd_constraint;
    const float dCDdX1y = -(f3*b1 + f3*b3 + f4*b2 + f4*b4)/cd_constraint;
    const float dCDdX2x = (f1*b1 + f2*b2)/cd_constraint; const float dCDdX2y = (f3*b1 + f4*b2)/cd_constraint;
    const float dCDdX3x = (f1*b3 + f2*b4)/cd_constraint; const float dCDdX3y = (f3*b3 + f4*b4)/cd_constraint;
    const float cd_grad_sum_sq =
      im1*(dCDdX1x*dCDdX1x + dCDdX1y*dCDdX1y) +
      im2*(dCDdX2x*dCDdX2x + dCDdX2y*dCDdX2y) +
      im3*(dCDdX3x*dCDdX3x + dCDdX3y*dCDdX3y);
    const float cd_pbd_lambda = -cd_constraint*cd_damping/cd_grad_sum_sq;
    const float dx1d = im1*cd_pbd_lambda*dCDdX1x; const float dy1d = im1*cd_pbd_lambda*dCDdX1y;
    const float dx2d = im2*cd_pbd_lambda*dCDdX2x; const float dy2d = im2*cd_pbd_lambda*dCDdX2y;
    const float dx3d = im3*cd_pbd_lambda*dCDdX3x; const float dy3d = im3*cd_pbd_lambda*dCDdX3y;

    pts[i1*dim+0] += dx1h + dx1d; pts[i1*dim+1] += dy1h + dy1d;
    pts[i2*dim+0] += dx2h + dx2d; pts[i2*dim+1] += dy2h + dy2d;
    pts[i3*dim+0] += dx3h + dx3d; pts[i3*dim+1] += dy3h + dy3d;
  }
}

void solve_length_constraints_tris_2d(float* pts, int num_pts, uint16_t* tris, int num_tris, float* orig,
    uint8_t* is_tri_active, float* inv_mass, float damping) {
  const int dim = 2;
  for (int i = 0; i < num_tris; i++) {
    if (is_tri_active[i] == 0) continue;
    const int i1 = tris[i*3+0]; const int i2 = tris[i*3+1]; const int i3 = tris[i*3+2];
    const vec2 op1 = (vec2){ orig[i1*dim+0], orig[i1*dim+1] };
    const vec2 op2 = (vec2){ orig[i2*dim+0], orig[i2*dim+1] };
    const vec2 op3 = (vec2){ orig[i3*dim+0], orig[i3*dim+1] };
    solve_length_constraints_of_pt_pairs_2d(pts, inv_mass, i1, i2, dist2(op1, op2), damping);
    solve_length_constraints_of_pt_pairs_2d(pts, inv_mass, i2, i3, dist2(op2, op3), damping);
    solve_length_constraints_of_pt_pairs_2d(pts, inv_mass, i3, i1, dist2(op3, op1), damping);
  }
}

void solve_length_constraints_tris_3d(float* pts, int num_pts, uint16_t* tris, int num_tris, float* orig,
    uint8_t* is_tri_active, float* inv_mass, float damping) {
  const int dim = 3;
  for (int i = 0; i < num_tris; i++) {
    if (is_tri_active[i] == 0) continue;
    const int i1 = tris[i*3+0]; const int i2 = tris[i*3+1]; const int i3 = tris[i*3+2];
    const vec3 op1 = (vec3){ orig[i1*dim+0], orig[i1*dim+1], orig[i1*dim+2] };
    const vec3 op2 = (vec3){ orig[i2*dim+0], orig[i2*dim+1], orig[i2*dim+2] };
    const vec3 op3 = (vec3){ orig[i3*dim+0], orig[i3*dim+1], orig[i3*dim+2] };
    solve_length_constraints_of_pt_pairs_3d(pts, inv_mass, i1, i2, dist3(op1, op2), damping);
    solve_length_constraints_of_pt_pairs_3d(pts, inv_mass, i2, i3, dist3(op2, op3), damping);
    solve_length_constraints_of_pt_pairs_3d(pts, inv_mass, i3, i1, dist3(op3, op1), damping);
  }
}

void solve_area_constraints_tris_2d(float* pts, int num_pts, uint16_t* tris, int num_tris, float* orig,
    uint8_t* is_tri_active, float* inv_mass, float damping) {
  const int dim = 2;
  for (int i = 0; i < num_tris; i++) {
    if (is_tri_active[i] == 0) continue;
    const int i1 = tris[i*3+0]; const int i2 = tris[i*3+1]; const int i3 = tris[i*3+2];
    const float  x1 = pts[i1*dim+0];  const float  x2 = pts[i2*dim+0];  const float  x3 = pts[i3*dim+0];
    const float  y1 = pts[i1*dim+1];  const float  y2 = pts[i2*dim+1];  const float  y3 = pts[i3*dim+1];
    const float ox1 = orig[i1*dim+0]; const float ox2 = orig[i2*dim+0]; const float ox3 = orig[i3*dim+0];
    const float oy1 = orig[i1*dim+1]; const float oy2 = orig[i2*dim+1]; const float oy3 = orig[i3*dim+1];
    const float im1 = inv_mass[i1];   const float im2 = inv_mass[i2];   const float im3 = inv_mass[i3];
    const float signed_area_orig = tri_area_2d_signed(ox1, oy1, ox2, oy2, ox3, oy3);
    const float signed_area_curr = tri_area_2d_signed( x1,  y1,  x2,  y2,  x3,  y3);
    const float coeff = .5 * (signed_area_curr < 0 ? -1. : 1.);
    const float dp1x = coeff*(y2-y3); const float dp1y = coeff*(x3-x2);
    const float dp2x = coeff*(y3-y1); const float dp2y = coeff*(x1-x3);
    const float dp3x = coeff*(y1-y2); const float dp3y = coeff*(x2-x1);
    const float grad_sum_sq = im1*(dp1x*dp1x + dp1y*dp1y) + im2*(dp2x*dp2x + dp2y*dp2y) + im3*(dp3x*dp3x + dp3y*dp3y);
    const float constraint = fabs(signed_area_curr) - fabs(signed_area_orig);
    const float lambda = -constraint*damping/grad_sum_sq;
    pts[i1*dim+0] += im1*lambda*dp1x; pts[i1*dim+1] += im1*lambda*dp1y;
    pts[i2*dim+0] += im2*lambda*dp2x; pts[i2*dim+1] += im2*lambda*dp2y;
    pts[i3*dim+0] += im3*lambda*dp3x; pts[i3*dim+1] += im3*lambda*dp3y;
  }
}

void solve_surface_area_constraints_tris_3d(float* pts, int num_pts, uint16_t* tris, int num_tris, float* orig,
  uint8_t* is_tri_active, float* inv_mass, float damping) {}

// could be optimized by not re-computing orig volume; and by not having to re-alloce 'grads' each time
void solve_volume_constraints_tris_3d(float* pts, int num_pts, uint16_t* tris, int num_tris, float* orig,
    uint8_t* is_tri_active, float* inv_mass, float damping, float pressure) {
  const int dim = 3;

  vec3 grads[num_pts];
  for (int i = 0; i < num_pts; i++) grads[i] = (vec3){ 0, 0, 0 };

  float curr_volume = 0.; float orig_volume = 0.;
  for (int i = 0; i < num_tris; i++) {
    //if (is_tri_active[i] == 0) return;
    const int i1 = tris[i*3+0]; const int i2 = tris[i*3+1]; const int i3 = tris[i*3+2];
    vec3 p1  = (vec3){ pts[i1*dim+0],  pts[i1*dim+1],  pts[i1*dim+2]  };
    vec3 p2  = (vec3){ pts[i2*dim+0],  pts[i2*dim+1],  pts[i2*dim+2]  };
    vec3 p3  = (vec3){ pts[i3*dim+0],  pts[i3*dim+1],  pts[i3*dim+2]  };
    vec3 op1 = (vec3){ orig[i1*dim+0], orig[i1*dim+1], orig[i1*dim+2] };
    vec3 op2 = (vec3){ orig[i2*dim+0], orig[i2*dim+1], orig[i2*dim+2] };
    vec3 op3 = (vec3){ orig[i3*dim+0], orig[i3*dim+1], orig[i3*dim+2] };
    curr_volume += 1./6. * dot3( p1, cross( p2,  p3));
    orig_volume += 1./6. * dot3(op1, cross(op2, op3));
    grads[i1] = add3(grads[i1], cross(p2, p3));
    grads[i2] = add3(grads[i2], cross(p3, p1));
    grads[i3] = add3(grads[i3], cross(p1, p2));
  }

  float grad_sum_sq = 0.;
  for (int i = 0; i < num_pts; i++)
    grad_sum_sq += inv_mass[i]*(grads[i].x*grads[i].x + grads[i].y*grads[i].y + grads[i].z*grads[i].z);

  const float constraint = curr_volume - orig_volume*pressure;
  const float lmd = -constraint*damping/grad_sum_sq;

  for (int i = 0; i < num_tris; i++) {
    const int i1 = tris[i*3+0]; const int i2 = tris[i*3+1]; const int i3 = tris[i*3+2];
    const float im1 = inv_mass[i1]; const float im2 = inv_mass[i2]; const float im3 = inv_mass[i3];
    pts[i1*dim+0] += im1*lmd*grads[i1].x; pts[i1*dim+1] += im1*lmd*grads[i1].y; pts[i1*dim+2] += im1*lmd*grads[i1].z;
    pts[i2*dim+0] += im2*lmd*grads[i2].x; pts[i2*dim+1] += im2*lmd*grads[i2].y; pts[i2*dim+2] += im2*lmd*grads[i2].z;
    pts[i3*dim+0] += im3*lmd*grads[i3].x; pts[i3*dim+1] += im3*lmd*grads[i3].y; pts[i3*dim+2] += im3*lmd*grads[i3].z;
  }
}

void apply_energy_based_breakage_tris_2d(float* pts, int num_pts, uint16_t* tris, int num_tris, float* orig,
    uint8_t* is_tri_active, float threshold1, float threshold2) {
  const int dim = 2;
  for (int i = 0; i < num_tris; i++) {
    if (is_tri_active[i] == 0) continue;
    const int i1 = tris[i*3+0]; const int i2 = tris[i*3+1]; const int i3 = tris[i*3+2];
    const float  x1 = pts[i1*dim+0];  const float  x2 = pts[i2*dim+0];  const float  x3 = pts[i3*dim+0];
    const float  y1 = pts[i1*dim+1];  const float  y2 = pts[i2*dim+1];  const float  y3 = pts[i3*dim+1];
    const float ox1 = orig[i1*dim+0]; const float ox2 = orig[i2*dim+0]; const float ox3 = orig[i3*dim+0];
    const float oy1 = orig[i1*dim+1]; const float oy2 = orig[i2*dim+1]; const float oy3 = orig[i3*dim+1];
    mat2x2 X    = (mat2x2){   x2-x1,   x3-x1,   y2-y1,   y3-y1 };
    mat2x2 Xhat = (mat2x2){ ox2-ox1, ox3-ox1, oy2-oy1, oy3-oy1 };
    const float detB = 1./det2(Xhat);
    const float b1 =  detB*Xhat.m22; const float b2 = -detB*Xhat.m12;
    const float b3 = -detB*Xhat.m21; const float b4 =  detB*Xhat.m11;
    const mat2x2 B = (mat2x2){ b1, b2, b3, b4 };
    const mat2x2 F = mat_mat_prod2x2(X, B);
    const float ch_energy = det2(F);
    const float cd_energy = sqrtf(F.m11*F.m11 + F.m12*F.m12 + F.m21*F.m21 + F.m22*F.m22);
    if (ch_energy > threshold1) is_tri_active[i] = 0;
    if (cd_energy > threshold2) is_tri_active[i] = 0;
  }
}

void apply_distance_based_breakage_tris_2d(float* pts, int num_pts, uint16_t* tris, int num_tris, float* orig,
  uint8_t* is_tri_active, float threshold) {}
void apply_distance_based_breakage_tris_3d(float* pts, int num_pts, uint16_t* tris, int num_tris, float* orig,
  uint8_t* is_tri_active, float threshold) {}

void test(float* arr) {
  for (int i = 0; i < 2; i++) arr[i] += 1;
}

