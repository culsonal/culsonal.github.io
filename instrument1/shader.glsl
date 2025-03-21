
uniform vec4 mouse;
uniform sampler2D img;
uniform int itr;
uniform float dt;

struct cell {
  float wall;
  float recording;
  float pressure;
  float pressure_prev;
  vec2 vel;
  vec2 vel_prev;
};

cell get_cell(ivec2 ixy) {
  vec4 t  = texelFetch(img, ixy, 0);
  vec2 wr = unpackHalf2x16(floatBitsToUint(t.x));
  vec2 p  = unpackHalf2x16(floatBitsToUint(t.y));
  vec2 vc = unpackHalf2x16(floatBitsToUint(t.z));
  vec2 vp = unpackHalf2x16(floatBitsToUint(t.w));
  return cell(wr.x, wr.y, p.x, p.y, vc, vp);
}

cell get_cell(vec2 xy) {
  vec4 t  = texelFetch(img, ivec2(xy*vec2(DIM)), 0);
  vec2 wr = unpackHalf2x16(floatBitsToUint(t.x));
  vec2 p  = unpackHalf2x16(floatBitsToUint(t.y));
  vec2 vc = unpackHalf2x16(floatBitsToUint(t.z));
  vec2 vp = unpackHalf2x16(floatBitsToUint(t.w));
  return cell(wr.x, wr.y, p.x, p.y, vc, vp);
}

vec4 cell_to_vec4(cell c) {
  float wr = uintBitsToFloat(packHalf2x16(vec2(c.wall, c.recording)));
  float p  = uintBitsToFloat(packHalf2x16(vec2(c.pressure, c.pressure_prev)));
  float vc = uintBitsToFloat(packHalf2x16(c.vel));
  float vp = uintBitsToFloat(packHalf2x16(c.vel_prev));
  return vec4(wr, p, vc, vp);
}

float sdf_box(vec2 p, vec2 b) {
  vec2 d = abs(p) - b;
  return length(max(d, vec2(0.))) + min(max(d.x, d.y), 0.);
}

float sdf_triangle(vec2 p) {
  const float k = sqrt(3.);
  p.x = abs(p.x) - 1.;
  p.y = p.y + 1. / k;
  if (p.x + k * p.y > 0.) p = vec2(p.x - k * p.y, -k * p.x - p.y) / 2.;
  p.x -= clamp(p.x, -2., 0.);
  return -length(p) * sign(p.y);
}

float sdf_capsule(vec2 p, vec2 a, vec2 b, float r) {
  vec2 pa = p - a;
  vec2 ba = b - a;
  float h = clamp(dot(pa, ba) / dot(ba, ba), 0., 1.);
  return length(pa - ba * h) - r;
}

float sdf_circle(vec2 p, float radius) {
  return length(p) - radius;
}

float sdf_subtract(float d1, float d2) { return max(d1, -d2); }
float sdf_union(float d1, float d2) { return min(d1, d2); }
float sdf_intersection(float d1, float d2) { return max(d1, d2); }

float flute_instrument(vec2 xy, vec2 dim) {
  float box_outside = sdf_box(xy*2.-1., dim);
  float box_inside = sdf_box(xy*2.-1., dim-vec2(.02));
  float hole = sdf_box(xy*2.-1.-vec2(dim.x, 0), vec2(dim.y, dim.y));
  float box = sdf_subtract(box_outside, box_inside);
  return sdf_subtract(box, hole);
}

float instrument1(vec2 xy) {
  float r = .4;
  return length(xy*2.-1.) > r && length(xy*2.-1.) < r+.029 && xy.x < .69 ? .9 : 0.;
}

float instrument2(vec2 xy) {
  return flute_instrument(xy, vec2(.6, .1)) < 0. ? 1. : 0.;
}

float instrument3(vec2 xy) {
  return flute_instrument(xy, vec2(.6, .3)) < 0. ? 1. : 0.;
}

float instrument4(vec2 xy) {
  return flute_instrument(xy, vec2(.6, .05)) < 0. ? 1. : 0.;
}

float instrument5(vec2 xy) {
  float tri_out = sdf_triangle((xy*2.-1.)*2.);
  float tri_in = sdf_triangle((xy*2.-1.)*2.1);
  float tri = sdf_subtract(tri_out, tri_in);
  return sdf_subtract(tri, sdf_box(xy*2.-1. + vec2(-.1, 0), vec2(.2, .2))) < 0. ? 1. : 0.;
}

float instrument6(vec2 xy) {
  float cap = sdf_capsule((xy*2.-1.), vec2(-.5, 0), vec2(.5, 0), .25);
  float c1 = sdf_circle((xy*2.-1.), .23);
  float c2 = sdf_circle((xy*2.-1.)+vec2( .44, 0), .23);
  float c3 = sdf_circle((xy*2.-1.)+vec2(-.44, 0), .23);
  float cap1 = sdf_subtract(cap, c1);
  float cap2 = sdf_subtract(cap1, c2);
  float cap3 = sdf_subtract(cap2, c3);
  float hole = sdf_circle(xy*2.-1.+vec2(-.67, 0), .12);
  float cap4 = sdf_subtract(cap3, hole);
  return cap4 < 0. ? 1. : 0.;
}

float instrument7(vec2 xy) {
  float flute = flute_instrument(xy, vec2(.6, .05));
  float hole1 = sdf_circle(xy*2.-1.+vec2(.5, -.05), .05);
  float hole2 = sdf_circle(xy*2.-1.+vec2(.3, -.05), .05);
  float f1 = sdf_subtract(flute, hole1);
  float f2 = sdf_subtract(f1,    hole2);
  return f2 < 0. ? 1. : 0.;
}

float instrument8(vec2 xy) {
  float flute = flute_instrument(xy, vec2(.6, .05));
  float hole1 = sdf_circle(xy*2.-1.+vec2(.5, -.05), .05);
  float hole2 = sdf_circle(xy*2.-1.+vec2(.3, -.05), .05);
  float hole3 = sdf_circle(xy*2.-1.+vec2(-.2, .05), .05);
  float f1 = sdf_subtract(flute, hole1);
  float f2 = sdf_subtract(f1,    hole2);
  float f3 = sdf_subtract(f2,    hole3);
  return f3 < 0. ? 1. : 0.;
}

float instrument9(vec2 xy) {
  return flute_instrument(xy, vec2(mouse.x, mouse.y)) < 0. ? 1. : 0.;
}

float instrument10(vec2 xy) {
  return flute_instrument(xy, vec2(.6, .03)) < 0. ? 1. : 0.;
}

cell init_state(ivec2 ixy, vec2 xy) {
  float wall = instrument8(xy); // todo try wall losses ie *.1 or whatever
  float init_p = 0.;
  float prev_p = init_p;
  vec2  init_v = vec2(0);
  vec2  prev_v = init_v;
  return cell(wall, 0., init_p, prev_p, init_v, prev_v);
}

float pml(vec2 xy) {
  float pml_thickness = .05;
  vec2 d = clamp(1.-min(xy, 1.-xy)/pml_thickness, 0., 1.);
  return 1.-dot(d, d);
}

cell update_state(ivec2 ixy, vec2 xy) {
  const float dt = .1;//.01;
  const float dx = 1.;
  cell c = get_cell(ixy);
  cell cxp = get_cell(ixy+ivec2( 1,  0));
  cell cxm = get_cell(ixy+ivec2(-1,  0));
  cell cyp = get_cell(ixy+ivec2( 0,  1));
  cell cym = get_cell(ixy+ivec2( 0, -1));
  float out_flow_x = (cxp.vel.x - c.vel.x)/(2.*dx);
  float out_flow_y = (cyp.vel.y - c.vel.y)/(2.*dx);
  float div = out_flow_x + out_flow_y;
  vec2 grad = vec2(
    (c.pressure - cxm.pressure)/(2.*dx),
    (c.pressure - cym.pressure)/(2.*dx)
  );

  float pressure_next = c.pressure_prev - dt*div;
  vec2 vel_next = c.vel_prev - dt*grad;
  //vel_next *= .9999;

  // input
  float d = distance(xy, vec2(0.5));
  float radius = 0.01;
  if (d < radius)
    pressure_next += .02*(1.-(d/radius))*sin(float(itr)/100. + 3.*sin(.0003*float(itr)));
    //pressure_next += .02*(1.-(d/radius))*sin(float(itr)/40.);

  // digested excitation into simpest form
  // (seems like adding in a pressure 'kick' isn't necessary after all)
  /*float d = distance(xy, vec2(.22, .5));
  float radius = 0.005;
  if (d < radius) {
    float p_mouth = .4;        // todo twiddle with this
    float p_bore = c.pressure;
    float delta_p = p_mouth - p_bore;
    if (delta_p > 0.) {
      vel_next += vec2(.03, .0) * delta_p;
      //pressure_next += .001 * delta_p;
    }
    //else vel_next = vec2(0);
  }*/

  // pml
  pressure_next *= pml(xy) * (1.-c.wall);
  vel_next *= pml(xy) * (1.-c.wall);

  // mic recording
  float recording_val;
  if (ixy.y*DIM.x + ixy.x == itr) recording_val = get_cell(mouse.xy).pressure;
  else recording_val = c.recording;

  return cell(
    c.wall, recording_val,
    pressure_next, c.pressure,
    vel_next, c.vel
  );
}

#ifdef RECORD
void record() {
  ivec2 ixy = ivec2(gl_FragCoord.xy);
  out_col = get_cell(ixy).recording;
}
#else
void update() {
  ivec2 ixy = ivec2(gl_FragCoord.xy);
  vec2  xy  = vec2(ixy)/vec2(DIM);
  if (itr == 0) out_col = cell_to_vec4(init_state(ixy, xy));
  else out_col = cell_to_vec4(update_state(ixy, xy));
  //else if (itr % 2 == 0) out_col = cell_to_vec4(update_state(ixy, xy));
  //else if (itr % 2 == 1) out_col = cell_to_vec4(advect_state(ixy, xy));
}

void draw() {
  ivec2 ixy = ivec2(gl_FragCoord.xy);
  vec2  xy  = vec2(ixy)/vec2(DIM);
  cell c = get_cell(ixy);
  if (distance(xy, mouse.xy) < .008) out_col = vec4(1);
  else if (c.wall > 0.) out_col = vec4(c.wall);
  else {
    float green = 0.;
    int rec_idx = int(xy.x*float(DIM.x*DIM.y));
    float mic_val_at_idx = get_cell(ivec2(rec_idx%DIM.x, rec_idx/DIM.x)).recording*.5+.5;
    mic_val_at_idx = clamp(mic_val_at_idx*.9, 0., .999);
    float graph_height = .3*float(DIM.y);
    float graph_y = mic_val_at_idx*graph_height;
    if (abs(graph_y-float(ixy.y)) < 1.) green = 1.-abs(graph_y-float(ixy.y));
   
    out_col = vec4(10.*vec3(max(0., -c.pressure), green, max(0., c.pressure)), 1);
    //out_col = vec4(1.*vec3(max(0., -c.pressure), c.recording, max(0., c.pressure)), 1);
    //out_col = vec4(vec3(max(0., -c.recording), length(c.vel), max(0., c.recording)), 1);
  }

  //else out_col = vec4(vec3(length(c.vel)), 1);
  //out_col = vec4(vec3(max(0., -c.recording), 0, max(0., c.recording)), 1);
}
#endif
