uniform mat4 matrix;

attribute vec2 v;
attribute vec2 tc;
varying vec2 v_texcoord;

void main() {
  gl_Position = matrix * vec4(v.x, v.y, 0.0, 1.0);
  v_texcoord = tc;
}
