#version 330

in vec2 v;
uniform vec4 colour;

out vec4 v_colour;

void main() {
  gl_Position = vec4(v.x, v.y, 0.0, 1.0);
  v_colour = colour;
}
