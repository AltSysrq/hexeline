uniform mat4 matrix;

attribute vec2 v;

void main() {
  gl_Position = matrix * vec4(v.x, v.y, 0.0, 1.0);
}
