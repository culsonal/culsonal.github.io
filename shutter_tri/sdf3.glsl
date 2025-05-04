
float sdf3_sphere(vec3 p, float r) {
  return length(p) - r;
}

float sdf3_torus(vec3 p, vec2 t) {
  vec2 q = vec2(length(p.xz) - t.x, p.y);
  return length(q) - t.y;
}

float sdf3_plane(vec3 p, vec4 n) {
  return dot(p, n.xyz) + n.w;
}

float sdf3_box(vec3 p, vec3 b) {
  vec3 q = abs(p) - b;
  return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float sdf3_octahedron(vec3 p, float s) {
  p = abs(p);
  float m = p.x + p.y + p.z - s;
  return m * (0.57735026919);
}

float sdf3_subtract(float d1, float d2) { return max(d1, -d2); }
float sdf3_union(float d1, float d2) { return min(d1, d2); }
float sdf3_intersection(float d1, float d2) { return max(d1, d2); }

