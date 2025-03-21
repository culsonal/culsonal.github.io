
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


float sdf_scene(vec2 xy) {
  float obj1 = sdf_subtract(
    sdf_box(xy, vec2(.95)),
    sdf_circle(xy-vec2(.9, 0), .5)
  );
  float obj2 = sdf_subtract(
    obj1,
    sdf_circle(xy+vec2(.9, 0), .5)
  );
  return obj2;
}



