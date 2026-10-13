// Camera.h -- a minimal mouse-driven orbit camera, built on plain GL1.1
// matrix-stack calls (glFrustum / gluLookAt-equivalent done by hand) so the
// renderer needs no GLU and no OpenGL function loader.
#pragma once
#include <cmath>

struct OrbitCamera {
    float yaw   = 0.9f;   // radians
    float pitch = 0.5f;
    float distance = 3.0f;
    float target[3] = {0,0,0};

    void applyProjection(float aspect) const {
        float fovy = 45.0f * 3.14159265f / 180.0f;
        float nearP = 0.05f, farP = 100.0f;
        float top = nearP * std::tan(fovy * 0.5f);
        float right = top * aspect;
        glFrustum(-right, right, -top, top, nearP, farP);
    }

    void applyView() const {
        float cx = target[0] + distance * std::cos(pitch) * std::sin(yaw);
        float cy = target[1] + distance * std::sin(pitch);
        float cz = target[2] + distance * std::cos(pitch) * std::cos(yaw);
        lookAt(cx, cy, cz, target[0], target[1], target[2], 0, 1, 0);
    }

    static void normalize3(float v[3]) {
        float len = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        if (len > 1e-9f) { v[0]/=len; v[1]/=len; v[2]/=len; }
    }
    static void cross3(const float a[3], const float b[3], float out[3]) {
        out[0] = a[1]*b[2]-a[2]*b[1];
        out[1] = a[2]*b[0]-a[0]*b[2];
        out[2] = a[0]*b[1]-a[1]*b[0];
    }

    static void lookAt(float ex,float ey,float ez, float cxp,float cyp,float czp, float ux,float uy,float uz) {
        float f[3] = { cxp-ex, cyp-ey, czp-ez };
        normalize3(f);
        float up[3] = { ux,uy,uz };
        float s[3]; cross3(f, up, s); normalize3(s);
        float u[3]; cross3(s, f, u);
        float m[16] = {
            s[0], u[0], -f[0], 0,
            s[1], u[1], -f[1], 0,
            s[2], u[2], -f[2], 0,
            0,    0,    0,     1
        };
        glMultMatrixf(m);
        glTranslatef(-ex, -ey, -ez);
    }
};
