// material class
#include "image_lib.h"
#ifndef HEADER_H
#define HEADER_H
#include <vector>
class Material {
public:
    float ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, ior, metal, rough;
    Material()
    {
        ar = ag = ab = dr = dg = db = sr = sg = sb = ns = tr = tg = tb = ior = metal = rough = 1;
    }
    Material(float a, float b, float c, float d, float e, float f, float g, float h, float i, float j, float k, float l, float m, float n,float o,float p)
    {
        ar = a;
        ag = b;
        ab = c;
        dr = d;
        dg = e;
        db = f;
        sr = g;
        sg = h;
        sb = i;
        ns = j;
        tr = k;
        tg = l;
        tb = m;
        ior = n;
        metal = o;
        rough = p;
    }

};
class Ray {
public:
    Point3D start;
    Line3D shoot;
    Dir3D dir;
    float intensity;
    Ray(Point3D start, Dir3D dir)
    {
        this->start = start;
        this->shoot = vee(start, dir).normalized();
        this->dir = dir.normalized();
        intensity = 1;

    }
};
class Hit {
public:
    Dir3D N;
    float intensity;
    float t; //distance from ray start
    Point3D interacts;
    Material material;
    int k; // index in scene
    int j; // 1 for sphere 2 for triangle
    Dir3D L;
    Dir3D V;
    Dir3D R;
    Hit()
    {
        intensity = -100;
    }
    Hit(Dir3D N, float intensity, Point3D interacts)
    {
        this->N = N.normalized();
        this->intensity = intensity;
        this->interacts = interacts;
        intensity = -100;
    }
    Hit(Dir3D N, Dir3D V,Dir3D R,Dir3D L,float intensity, Point3D interacts)
    {
        this->N = N.normalized();
        this->intensity = intensity;
        this->interacts = interacts;
        this->N = N;
        this->V = V;
        this->L = L;
        intensity = -100;
    }
};
class PointLights {
public:
    Color light;
    Point3D center;
    PointLights(Color color, Point3D position)
    {
        light = color;
        center = position;
    }
};
class DirectionLights {
public:
    Color light;
    Dir3D direction;
    DirectionLights(Color color, Dir3D dir)
    {
        light = color;
        direction = dir;
    }
};

class Triangles {
public:
    Point3D v1;
    Point3D v2;
    Point3D v3;
    Dir3D n1;
    Dir3D n2;
    Dir3D n3;
    Material material;
    Triangles(Point3D v1, Point3D v2, Point3D v3, Dir3D n1)
    {
        this->v1 = v1;
        this->v2 = v2;
        this->v3 = v3;
        this->n1 = n1;
        this->n2 = n1;
        this->n3 = n1;
    }
    Triangles(Point3D v1,Point3D v2,Point3D v3,Dir3D n1,Dir3D n2,Dir3D n3)
    {
        this->v1 = v1;
        this->v2 = v2;
        this->v3 = v3;
        this->n1 = n1;
        this->n2 = n2;
        this->n3 = n3;
    }
};
class Boxes
{

};
class Planes
{

};
class Objects
{
public:
    Point3D pos;
    Material material;

};
class Spheres{
public:
    float radius;
    Point3D pos;
    Material material;
    Spheres()
    {
        pos = Point3D(0, 0, 2);
        radius = 1;
    }
    Spheres(Point3D center, float radius, Material mat)
    {
        pos = center;
        this->radius = radius;
        material = mat;
    }
    

};
class Scene
{
public:
    std::vector<Spheres> spheres;
    std::vector<Triangles> triangles;
    std::vector<Boxes> boxes;
    std::vector<Planes> planes;
    std::vector<Objects> objects;
    std::vector<PointLights> pointLights;
    std::vector<DirectionLights> directLights;
};

#endif