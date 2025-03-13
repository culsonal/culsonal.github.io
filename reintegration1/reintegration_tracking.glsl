uniform sampler2D tex;
uniform vec4 mouse;
uniform int itr;
out vec4 out_col;

struct particle {
  float mass;
  vec2 vel;
  vec2 cm;
};

particle unpack_particle(ivec2 xy) {
  vec4 t = texelFetch(tex, xy, 0);
  return particle(
    t.x,
    unpackSnorm2x16(floatBitsToUint(t.y)),
    unpackUnorm2x16(floatBitsToUint(t.z))
  );
}

vec4 pack_particle(particle p) {
  return vec4(
    p.mass,
    uintBitsToFloat(packSnorm2x16(p.vel)),
    uintBitsToFloat(packUnorm2x16(p.cm)),
    1
  );
}

vec3 area_overlap(vec2 dxy, float K) {
  vec2 omin = clamp(dxy-K*.5, 0., 1.);
  vec2 omax = clamp(dxy+K*.5, 0., 1.);
  return vec3(.5*(omin + omax), (omax.x - omin.x)*(omax.y - omin.y)/(K*K));
}

particle init_state(ivec2 ixy) {
  //float(ixy.x%3 == 0 && ixy.y%3 == 1) * hash12(vec2(ixy)) * 5.,
  //ixy.y < DIM.y-20 ? 0. : hash12(vec2(ixy))*2.,
  return particle(
    length(vec2(ixy)/float(DIM)-vec2(.5, .35)) > .15 ? 0. : hash12(vec2(ixy))*2.,
    hash22(vec2(ixy))*4.-2., //vec2(0, 1),
    vec2(.5)
  );
}

particle advect_particle(ivec2 ixy, float dt) {
  particle p_next = particle(0., vec2(0), vec2(0));

  for (int dx = -2; dx <= 2; dx++) for (int dy = -2; dy <= 2; dy++) {
    vec2 dxy = vec2(dx, dy);
    particle np = unpack_particle(ixy + ivec2(dxy));
    vec2 neigh_next_pos = dxy + np.cm + dt*np.vel;

    float diffusion = 1.25;
    vec3 A = area_overlap(neigh_next_pos, diffusion);
    vec2  overlap_cm   = A.xy;
    float overlap_mass = A.z * np.mass;
    p_next.mass += overlap_mass;
    p_next.vel  += np.vel*overlap_mass;
    p_next.cm   += overlap_cm*overlap_mass;
  }

  if (p_next.mass > 0.) {
    p_next.vel /= p_next.mass;
    p_next.cm  /= p_next.mass;
  }

  return p_next;
}

void update_via_sph_pressure(inout particle p, ivec2 ixy, float dt) {
  vec2 acceleration = vec2(0, 0);
  for (int dx = -2; dx <= 2; dx++) for (int dy = -2; dy <= 2; dy++) {
    vec2 dxy = vec2(dx, dy);
    particle np = unpack_particle(ixy + ivec2(dxy));

    const float k = .01;
    const float volume = 1.;
    float rho_0 =  .8;
    float rho_i =  p.mass/volume;
    float rho_j = np.mass/volume;
    float P_i = k*rho_i*(rho_i-rho_0);
    float P_j = k*rho_j*(rho_j-rho_0);
    vec2 diff = dxy + np.cm  -  p.cm;
    float sm_kern = 2.*exp(-dot(diff, diff));
    if (rho_i > 0. && rho_j > 0.)
      acceleration -= np.mass*(P_i/(rho_i*rho_i) + P_j/(rho_j*rho_j))*sm_kern*diff;
  }

  p.vel += dt*acceleration;
}

particle update_particle_physics(ivec2 ixy, float dt) {
  particle p = unpack_particle(ixy);
  update_via_sph_pressure(p, ixy, dt);

  //p.vel.y -= .001;
  if (ixy.x == 1       && p.vel.x < 0.) p.vel.x *= -1.;
  if (ixy.y == 1       && p.vel.y < 0.) p.vel.y *= -1.;
  if (ixy.x == DIM.x-1 && p.vel.x > 0.) p.vel.x *= -1.;
  if (ixy.y == DIM.y-1 && p.vel.y > 0.) p.vel.y *= -1.;
  if (distance(vec2(ixy)/vec2(DIM), mouse.xy) < .05) p.vel += mouse.zw*10.;
  if (length(p.vel*dt) > 1.) p.vel = p.vel/length(p.vel*dt);
  return p;
}

void update() {
  const float dt = .5;
  ivec2 ixy = ivec2(gl_FragCoord.xy);
  if (itr == 0)          out_col = pack_particle(init_state(ixy));
  else if (itr % 2 == 0) out_col = pack_particle(advect_particle(ixy, dt));
  else                   out_col = pack_particle(update_particle_physics(ixy, dt));
}

void draw() {
  /*ivec2 ixy = ivec2(gl_FragCoord.xy);
  vec2  uv  = gl_FragCoord.xy/float(DIM);
  float mass = texelFetch(tex, ixy, 0).x*.01;
  vec3 dfdx = dFdx(vec3(uv, mass));
  vec3 dfdy = dFdy(vec3(uv, mass));
  vec3 normal = normalize(cross(dfdx, dfdy));
  vec3 light_pos = vec3(.5, .5, .4);
  vec3 light_dir = normalize(light_pos - vec3(uv, 0));
  float shade = max(0., dot(normal, light_dir));
  out_col = vec4(vec3(shade*.8), 1);*/

  ivec2 ixy = ivec2(gl_FragCoord.xy);
  particle p = unpack_particle(ixy);
  //out_col = vec4(p.mass, (p.vel*.5+.5)*p.mass, 1);
  out_col = vec4(vec3(p.mass), 1);
}

