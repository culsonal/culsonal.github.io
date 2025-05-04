
uniform vec3 mouse;
uniform float itr;
uniform float focal_length;
out vec4 out_col;

float T = 0.;

float sdf_scene(vec3 p) {
  //vec3 sc = vec3(0, 0, .5); p += sc;
  vec3 sc = vec3(sin(T*6.28)*.15, cos(T*6.28)*.15, .5); p += sc;
  float s1 = sdf3_sphere(p+vec3(.05,  .1, .2), .1);
  float s2 = sdf3_sphere(p+vec3(.01, .08, .31), .04);
  float s3 = sdf3_subtract(s1, s2);
  //float s4 = sdf3_plane(p-vec3(0, 1, 0), vec4(0, -1, 0, 0));
  //float s5 = sdf3_union(s3, s4);
  //return s4;
  return s3;
}

vec3 sdf_normal(vec3 p) {
  float eps = 0.0001;
  return normalize(vec3(
    sdf_scene(p+vec3(eps, 0, 0)) - sdf_scene(p-vec3(eps, 0, 0)),
    sdf_scene(p+vec3(0, eps, 0)) - sdf_scene(p-vec3(0, eps, 0)),
    sdf_scene(p+vec3(0, 0, eps)) - sdf_scene(p-vec3(0, 0, eps))
  ));
}

int shutter_mask(vec2 xy, int itr, int max_itr) {
  float t = float(itr)/float(max_itr);
  float size = exp((1.-t)*5.);
  //float s = sdf2_box(xy*size, vec2(1., 1.));
  float s = sdf2_circle(xy*size, 1.5);
  //float s = sdf2_triangle(xy*size*.4);
  //float s = sdf2_capsule(xy*size+vec2(1, 0), vec2(2., 0), vec2(0), 1.);
  // todo some union/intersection geometries
  return s > 0. ? 0 : 1;
}

float ray_march(vec3 ro, vec3 rd) {
  float t = 0.;
  for (int i = 0; i < 500; i++) {
    float d = sdf_scene(ro + t*rd);
    if (d < .00001) return t;
    t += d;
    if (t > 50.) break;
  }
  return -1.;
}

float shading(vec3 p, vec3 normal) {
  //vec3 light = vec3(0. + sin(itr/50.), 0, -2.);
  vec3 light = vec3(.8, 0, -2.);
  vec3 L = normalize(light - p);
  return max(dot(normal, L), 0.);
}

vec3 compute_color(vec2 xy) {
  vec3 ray_origin = vec3(xy, -focal_length);
  vec3 pinhole = vec3(0);
  vec3 ray_dir = normalize(pinhole - ray_origin);
  float t = ray_march(ray_origin, ray_dir);
  if (t < 0.) return vec3(0);
  vec3 p = ray_origin + t*ray_dir;
  vec3 n = sdf_normal(p);
  return vec3(1) * shading(p, n);

  //return vec3(1.-distance(ray_origin, p)/10.);
  //return vec3(1.-distance(ray_origin, p)/2.);
  //return vec3(distance(ray_origin, p)/2.);
}

vec3 compute_color_over_shutter_seq(vec2 xy) {
  int shutter_speed = 100;
  vec3 tot_col = vec3(0, 0, 0);
  float num_times_active = 0.;
  for (int i = 0; i < shutter_speed; i++) {
    int mask = shutter_mask(xy, i, shutter_speed);
    if (mask == 0) continue;
    num_times_active += 1.;
    T = (itr+float(i))/float(shutter_speed) * .2; // animation speed control
    tot_col += compute_color(xy);
  }

  //float exposure = 1./float(shutter_speed) * 5.;
  //return tot_col * exposure;
  return tot_col / num_times_active;
}

void draw() {
  //T = float(itr/100.);

  vec2 xy = vec2(gl_FragCoord.xy)/vec2(DIM)*2.-1.;
  out_col = vec4(compute_color_over_shutter_seq(xy), 1);

  //int mask = shutter_mask(xy, itr, 100.);
  //int mask = shutter_mask(xy, 98., 100.);
  //if (mask == 0) out_col = vec4(vec3(.1), 1);
  //if (mask == 1) out_col = vec4(compute_color(xy), 1);
}


