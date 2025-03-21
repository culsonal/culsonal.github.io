uniform vec3 mouse;
uniform vec2 mouse_vel;
uniform int num_sim_steps;
uniform sampler2D img;
uniform int itr;

struct cell {
  float is_active;
  float recording;
  vec2 off;
  vec2 vel;
};

cell get_cell(ivec2 ixy) {
  vec4 t = texelFetch(img, ixy, 0);
  vec2 ar  = unpackHalf2x16(floatBitsToUint(t.x));
  vec2 off = unpackHalf2x16(floatBitsToUint(t.y));
  vec2 vel = unpackHalf2x16(floatBitsToUint(t.z));
  return cell(ar.x, ar.y, off, vel);
}

cell get_cell(vec2 xy) {
  return get_cell(ivec2(xy*vec2(DIM)));
}

vec4 cell_to_vec4(cell c) {
  return vec4(
    uintBitsToFloat(packHalf2x16(vec2(c.is_active, c.recording))),
    uintBitsToFloat(packHalf2x16(c.off)),
    uintBitsToFloat(packHalf2x16(c.vel)),
    1.
  );
}

cell init_state(ivec2 ixy, vec2 xy) {
  return cell(sdf_scene(xy*2.-1.) < 0. ? 1. : 0., 0.,
    vec2(0),
    vec2(0));
    //length(xy*2.-1.) < .08 ? vec2(-.08, 0) : vec2(0));
}

vec2 edge_length_force(ivec2 ixy) {
  cell c = get_cell(ixy);
  if (c.is_active == 0.) return vec2(0);
  cell cxm = get_cell(ixy+ivec2(-1,  0));
  cell cxp = get_cell(ixy+ivec2( 1,  0));
  cell cym = get_cell(ixy+ivec2( 0, -1));
  cell cyp = get_cell(ixy+ivec2( 0,  1));
  float spring_k = .5;
  vec2 fxm = spring_k*(cxm.off - c.off);
  vec2 fxp = spring_k*(cxp.off - c.off);
  vec2 fym = spring_k*(cym.off - c.off);
  vec2 fyp = spring_k*(cyp.off - c.off);
  return fxm + fxp + fym + fyp;
}

vec2 mouse_vel_if_active(vec2 xy) {
  if (mouse.z > .5 && distance(xy, mouse.xy) < .08) return mouse_vel*.1;
  else return vec2(0);
}

float capture_recording(ivec2 ixy, float curr_rec) {
  //if (ixy.y*DIM.x + ixy.x == itr) return length(get_cell(DIM/2).off);
  //if (ixy.y*DIM.x + ixy.x == itr) return get_cell(DIM/2).vel.y;
  if (ixy.y*DIM.x + ixy.x == itr) return get_cell(DIM/2).off.x;
  else return curr_rec;
}

cell update_state(ivec2 ixy, vec2 xy) {
  const float vel_damping = .999;
  const float dt = .9;
  const float m = 1.;
  cell c = get_cell(ixy);
  vec2 f = edge_length_force(ixy);
  vec2 next_vel = (c.vel + dt*f/m + mouse_vel_if_active(xy)) * vel_damping;
  vec2 next_off =  c.off + dt*next_vel;
  float rec = capture_recording(ixy, c.recording);
  return cell(c.is_active, rec, next_off, next_vel);
}

#ifdef RECORD
void record() {
  ivec2 ixy = ivec2(gl_FragCoord.xy);
  out_col = get_cell(ixy).recording;
}
#else
void update() {
  ivec2 ixy = ivec2(gl_FragCoord.xy);
  vec2  xy = gl_FragCoord.xy/vec2(DIM);
  if (itr == 0) out_col = cell_to_vec4(init_state(ixy, xy));
  else          out_col = cell_to_vec4(update_state(ixy, xy));
}

void draw() {
  ivec2 ixy = ivec2(gl_FragCoord.xy);
  vec2  xy = gl_FragCoord.xy/vec2(DIM);
  cell c = get_cell(ixy);
  bool cxm = get_cell(ixy+ivec2(-1,  0)).is_active == 0.;
  bool cxp = get_cell(ixy+ivec2( 1,  0)).is_active == 0.;
  bool cym = get_cell(ixy+ivec2( 0, -1)).is_active == 0.;
  bool cyp = get_cell(ixy+ivec2( 0,  1)).is_active == 0.;
  if (distance(xy, mouse.xy) < .004) out_col = vec4(1);
  else if (c.is_active == 1. && (cxm || cxp || cym || cyp)) out_col = vec4(1);
  else {
    //int rec_idx = int(xy.x*float(DIM.x*DIM.y));
    int rec_idx = int(xy.x*float(num_sim_steps));
    float rec = get_cell(ivec2(rec_idx%DIM.x, rec_idx/DIM.x)).recording*.5+.5;
    rec = clamp(rec*.9, 0., .999);
    float graph_height = .3*float(DIM.y);
    float graph_y = rec*graph_height;
    float graph_col = 0.;
    if (abs(graph_y-float(ixy.y)) < 1.) graph_col = 1.-abs(graph_y-float(ixy.y));

    float o = c.is_active == 0. ? 0. : length(c.off);
    out_col = vec4(o, graph_col, o, 1);
  }
}
#endif
