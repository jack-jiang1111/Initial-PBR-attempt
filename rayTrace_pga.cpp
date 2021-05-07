#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//Images Lib includes:
#define STB_IMAGE_IMPLEMENTATION //only place once in one .cpp file
#define STB_IMAGE_WRITE_IMPLEMENTATION //only place once in one .cpp files
#define M_PI 3.1415926
#include "image_lib.h" //Defines an image class and a color class

//#3D PGA
#include "PGA_3D.h"

//High resolution timer
#include <chrono>

//Scene file parser
#include "parse_pga.h"
#include <vector>
Color EvaluateRayTree(Scene scene, Ray ray, int depth);
bool FindIntersection(Scene scene, Ray ray, Hit& hit,int k, int j);
Color ApplyLightingModel(Scene scene, Ray ray, Hit hit, int depth);
Color diffuseContribution(PointLights L, Hit hit);
Color SpecularContribtion(PointLights L, Hit hit);
Color diffuseContribution(DirectionLights L, Hit hit);
Color SpecularContribtion(DirectionLights L, Hit hit);
Color ReflectAndRefractContribtion(Color reflect, Color refract, Hit hit);
Color ReflectContribtion(Color reflect, Hit hit);

Ray Reflect(Ray ray, Dir3D norm, Point3D start);
Ray Refract(Ray ray, Hit hit);
Color Ambient(Color ambient, Material mat);
std::vector<Ray> ReflectWithParameter(Ray input, Hit hit);


float raySphereIntersect(Ray ray, Spheres sphere, Hit& hit) {
    Line3D rayLine = ray.shoot;
    Point3D sphereCenter = sphere.pos;
    float sphereRadius = sphere.radius;
    Point3D rayStart = ray.start;
    Point3D projPoint = dot(rayLine, sphereCenter) * rayLine;      //Project to find closest point between circle center and line [proj(sphereCenter,rayLine);]
    float distSqr = projPoint.distToSqr(sphereCenter);          //Point-line distance (squared)
    float d2 = distSqr / (sphereRadius * sphereRadius);             //If distance is larger than radius, then...
    if (d2 > 1)
    {
        //std::cout << "!!!!!!!" << std::endl;
        return -1;                                   //... the ray missed the sphere
    }
    float w = sphereRadius * sqrt(1 - d2);                          //Pythagorean theorem to determine dist between proj point and intersection points
    Point3D p1 = projPoint + rayLine.dir() * w;                   //Add/subtract above distance to find hit points
    Point3D p2 = projPoint - rayLine.dir() * w;

    if (dot((p1 - rayStart), rayLine.dir()) >= 0) //Is the first point in same direction as the ray line?
    {
        hit.material = sphere.material;
        hit.N = (p1 - sphere.pos).normalized();
        hit.interacts = p1;
        hit.t = p1.distTo(ray.start);
        return hit.t;     
    }
    if (dot((p2 - rayStart), rayLine.dir()) >= 0)  //Is the second point in same direction as the ray line?
    {
        hit.material = sphere.material;
        hit.N = (p2 - sphere.pos).normalized();
        hit.interacts = p2;
        hit.t = p2.distTo(ray.start);
        return hit.t;     
    }
    //std::cout << "??" << std::endl;
    return -1;
}
bool sameside(Point3D p1, Point3D p2, Point3D a, Point3D b)
{
    Dir3D cp1 = cross(b - a, p1 - a);
    Dir3D cp2 = cross(b - a, p2 - a);
    return dot(cp1, cp2) >= 0;
}
bool inside_triangle(Point3D interact_point, Triangles triangle)
{
    return sameside(interact_point, triangle.v1, triangle.v2, triangle.v3) && sameside(interact_point, triangle.v2, triangle.v1, triangle.v3) && sameside(interact_point, triangle.v3, triangle.v1, triangle.v2);
}
float rayTriangleIntersect(Ray ray, Triangles triangle, Hit &hit) { // if ray and sphere intersect
    Dir3D dir = triangle.v1 - ray.start;
    float l_dot_n = dot(ray.dir, triangle.n1);
    if (l_dot_n == 0) // parallal to triangle
    {
        return -1;
    }

    float dis = dot(dir,triangle.n1)/l_dot_n; // this is the distance between ray start and triangle plane

    Point3D interact_point = ray.start + dis * ray.dir;
    if (inside_triangle(interact_point, triangle))
    {
        hit.interacts = interact_point;
        hit.N = triangle.n1;
        hit.material = triangle.material;
        hit.t = dis;
        return dis;
    }
    else
    {
        return - 1;
    }

    // method 2:
    /*
    Plane3D plane = wedge(triangle.v1, triangle.v2, triangle.v3);
    float d = fabs(plane.w) / plane.magnitudeSqr;
    float t = -1 * (dot(ray.start, triangle.n1) + d) / dot(ray.dir, triangle.n1);
    Point3D interact = ray.start + t * ray.dir;
    */
}

int main(int argc, char** argv){
  
  //Read command line paramaters to get scene file
  if (argc != 2){
     std::cout << "Usage: ./a.out scenefile\n";
     return(0);
  }
  std::string secenFileName = argv[1];

  //Parse Scene File
  parseSceneFile(secenFileName);

  float imgW = img_width, imgH = img_height;
  float halfW = imgW/2, halfH = imgH/2;
  float d = halfH / tanf(halfAngleVFOV * (M_PI / 180.0f));
  // d = halfW / tanf(halfAngleVFOV * (M_PI / 180.0f));
  Color back = backgroundcolor;
  Image outputImg = Image(img_width,img_height);
  auto t_start = std::chrono::high_resolution_clock::now();


  for (int i = 0; i < img_width; i++){
    for (int j = 0; j < img_height; j++){
      //TODO: In what way does this assumes the basis is orthonormal?
        float u,v;
        if (sampling == 0) // jittered supersampling
        {
            float rand1 = (float)rand() / (RAND_MAX);
            float rand2 = (float)rand() / (RAND_MAX);
            u = (halfW - (imgW) * ((i + rand1) / imgW));
            v = (halfH - (imgH) * ((j + rand2) / imgH));
        }
        else if (sampling == 1) //adaptive supersampling
        {

        }
        else// basic sampling
        {
            u = (halfW - (imgW) * ((i + 0.5) / imgW));
            v = (halfH - (imgH) * ((j + 0.5) / imgH));
        }
      
      Point3D p = eye - d*forward + u*right + v*up;
      Dir3D rayDir = p - eye;
      Ray view(eye, rayDir);
      Color color = EvaluateRayTree(scene, view, max_depth);
      //color.print();
      outputImg.setPixel(i, j, color);
      //outputImg.setPixel(i,j, Color(fabs(i/imgW),fabs(j/imgH),fabs(0))); //TODO: Try this, what is it visualizing?
    }
  }
  auto t_end = std::chrono::high_resolution_clock::now();
  printf("Rendering took %.2f ms\n",std::chrono::duration<double, std::milli>(t_end-t_start).count());

  outputImg.write(imgName.c_str());
  return 0;
}

Color EvaluateRayTree(Scene scene, Ray ray,int depth) {
    bool hit_something;
    Hit HitInformation; // structure containing hit point, normal, etc
    hit_something = FindIntersection(scene, ray, HitInformation,-1,-1);
    
    if (hit_something) {
        
        return ApplyLightingModel(scene,ray, HitInformation ,depth);
    }
    else {
        
        return backgroundcolor;
    }
}
bool FindIntersection(Scene scene, Ray ray, Hit &hit, int k, int j){ 
    // int k is the index in scene array
    // int j  == 1 -> sphere, int j == 2 ->triangle
    float min_t = INT_MAX;
    bool hits = false;
    Hit defaulthit;

    for (int i = 0; i < scene.spheres.size(); i++)
    {
        if (j == 1 && i == k)
        {
            continue; // do not check itself
        }
        Spheres hitsphere = scene.spheres[i];
        float d = raySphereIntersect(ray, hitsphere, defaulthit);

        if (d > 0 && d < min_t)
        {
            min_t = d;
            hits = true;
            hit = defaulthit;
            hit.intensity = 1/d;
            hit.j = 1;
            hit.k = i;
        }
    }
    for (int i = 0; i < scene.triangles.size(); i++)// check 
    {
        if (j == 2 && i == k)
        {
            continue; // do not check itself
        }
        Triangles hitTriangle = scene.triangles[i];
        float d = rayTriangleIntersect(ray, hitTriangle,  defaulthit);
        if (d > 0 && d < min_t)
        {
            min_t = d;
            hits = true;
            hit = defaulthit;
            hit.intensity = 1/d;
            hit.j = 2;
            hit.k = i;
        }
    }
    //  mstd::cout << hit.intensity << std::endl;
    return hits;
}
Color ApplyLightingModel(Scene scene, Ray ray, Hit hit,int depth)
{
    Color contribution = backgroundcolor;
    for (int i = 0; i < scene.pointLights.size(); i++) // sum over all point light
    {
        Dir3D dir = (scene.pointLights[i].center - hit.interacts).normalized();
        Ray shadow(hit.interacts + dir * 0.01, dir);
        Hit shadow_hit;
        bool blocked = FindIntersection(scene, shadow, shadow_hit,hit.k,hit.j);

        Dir3D d = scene.pointLights[i].center - hit.interacts;
        
        if (blocked && shadow_hit.t<d.magnitude()) {
            //contribution = contribution + Color(1, 1, 1);
            continue; // we¡¯re in shadow, on to the next light;
        }
        else
        {
            //std::cout << contribution.r << " " << contribution.g << " " << contribution.b << std::endl;
            float percent = 1; // how different percent contributing the color
            contribution = contribution + diffuseContribution(scene.pointLights[i], hit)*percent;
            contribution = contribution + SpecularContribtion(scene.pointLights[i], hit)* (1 - percent);
            //std::cout << contribution.r << " " << contribution.g << " " << contribution.b << std::endl << std::endl;
        }
        
    }
    for (int i = 0; i < scene.directLights.size(); i++) //sum over direction light
    {
        Dir3D dir = scene.directLights[i].direction*(-1);
        // normal 
        Ray shadow(hit.interacts+hit.N*0.01, dir);
        Hit shadow_hit;
        //std::cout << shadow_hit.intensity << std::endl;
        bool blocked = FindIntersection(scene, shadow, shadow_hit, hit.k, hit.j);
        //std::cout << shadow_hit.intensity << " " << blocked << std::endl << std::endl; 
        if (!blocked) {
            
            //std::cout << contribution.r << " " << contribution.g << " " << contribution.b << std::endl;
            contribution = contribution + diffuseContribution(scene.directLights[i], hit);
            contribution = contribution + SpecularContribtion(scene.directLights[i], hit);
            //std::cout << contribution.r << " " << contribution.g << " " << contribution.b << std::endl << std::endl;
        }
        else
        {
            //std::cout << shadow_hit.intensity << " " << blocked << std::endl << std::endl;
            contribution = contribution + Color(1, 1, 1);
        }
    }
    if (depth != 0)
    {
        // energy conserve case
        /*
        std::vector<Ray> ray_queue = ReflectWithParameter(ray, hit);
        for (int i = 0; i < ray_queue.size(); i++)
        {
            contribution = contribution + ray_queque[i].intensity*ReflectContribtion(EvaluateRayTree(scene, ray_queue[i], depth - 1), hit);
        }
        */
        Ray mirror = Reflect(ray, hit.N, hit.interacts);
        Color reflect = EvaluateRayTree(scene, mirror, depth - 1);
        Ray glass = Refract(ray, hit);
        Color refract = EvaluateRayTree(scene, glass, depth - 1);
        contribution = contribution + ReflectAndRefractContribtion(reflect, refract, hit);
    }
    contribution = contribution + Ambient(ambient,hit.material); 
    return contribution;
}
Color Ambient(Color ambient, Material mat)
{
    return Color(ambient.r * mat.ar, ambient.g * mat.ag, ambient.b * mat.ab);
}
Color diffuseContribution(PointLights L, Hit hit)
{
    Dir3D light = (L.center - hit.interacts);
    float intensity = 1 / (0.1 * light.magnitude()); 
    light = light.normalized();
    float factor = dot(hit.N,light);
    
    //std::cout << factor << std::endl;
    if (factor < 0)
    {
        factor = 0;
    }
    return Color(hit.material.dr * factor* intensity, hit.material.dg * factor* intensity, hit.material.db * factor* intensity);
}
Color diffuseContribution(DirectionLights L, Hit hit)
{
    Dir3D light = L.direction*(-1);
    float intensity = 0.75;
    float factor = dot(hit.N, light);

    //std::cout << factor << std::endl;
    if (factor < 0)
    {
        factor = 0;
    }
    return Color(hit.material.dr * factor * intensity, hit.material.dg * factor * intensity, hit.material.db * factor * intensity);
}
Color SpecularContribtion(PointLights L, Hit hit)
{
    Dir3D V = (eye - hit.interacts).normalized();
    Dir3D light = ( hit.interacts- L.center); // L 
    float intensity = 1 / (0.1 * light.magnitude()); // intensity
    light = light.normalized(); // normalize
    Dir3D R = (light - 2 * dot(light, hit.N) * hit.N).normalized();
    float factor = pow(dot(V, R), hit.material.ns);
    return Color(hit.material.sr * factor * intensity, hit.material.sg * factor * intensity, hit.material.sb * factor * intensity);
}
Color SpecularContribtion(DirectionLights L, Hit hit)
{
    Dir3D V = (eye - hit.interacts).normalized();
    Dir3D light = L.direction ;
    float intensity = 2; // intensity
    Dir3D R = (light - 2 * dot(light, hit.N) * hit.N).normalized();
    float factor = pow(dot(V, R), hit.material.ns);
    return Color(hit.material.sr * factor * intensity, hit.material.sg * factor * intensity, hit.material.sb * factor * intensity);
}
Color ReflectContribtion(Color reflect, Hit hit)
{
    float red = reflect.r * hit.material.sr*hit.intensity;
    float green = reflect.g * hit.material.sg*hit.intensity;
    float blue = reflect.b * hit.material.sb*hit.intensity;
    return Color(red, green, blue);
}
Color ReflectAndRefractContribtion(Color reflect,Color refract, Hit hit)
{
    float red = reflect.r * hit.material.sr + refract.r * hit.material.tr;
    float green = reflect.g * hit.material.sg + refract.g * hit.material.tg;
    float blue = reflect.b * hit.material.sb + refract.b * hit.material.tb;
    return Color(red, green, blue);
}
Ray Reflect(Ray ray, Dir3D norm, Point3D start)
{
    Dir3D R = ray.dir - 2 * dot(ray.dir, norm) * norm;
    R = R.normalized();
    return Ray(start+0.01*R, R);
}
Ray Refract(Ray ray, Hit hit)
{
    float in_or_out = dot(hit.N, ray.dir);
    float ni, nr;
    float anglei, angler;
    nr = hit.material.ior;
    ni = 1;
    anglei = acos(dot(hit.N, ray.dir) / ray.dir.magnitude());
    angler = asin(ni * sin(anglei) / nr);
    /*
    if (in_or_out > 0) // out
    {
        ni = hit.material.ior;
        nr = 1;
        anglei
    }
    else // in
    {
        nr = hit.material.ior;
        ni = 1;
        anglei = acos(dot(hit.N, ray.dir) / ray.dir.magnitude());
        angler = asin(ni * sin(anglei) / nr);
    }*/
    // orthonormal case?

    Dir3D T = (ni / nr * cos(anglei) - cos(angler)) * hit.N - ni / nr * ray.dir;
    return Ray(hit.interacts, T.normalized());
}

std::vector<Ray> ReflectWithParameter(Ray input, Hit hit)
{
    // take an input ray, based on parameter( roughness of material), return mutiple ray, give them proper intensity
    std::vector<Ray> result;
    // the reflection intensity is equal to metal parameter
    // a higher metal value means reflecting more light
    float reflection_intensity = hit.material.metal;
    // 1- reflection_intensity is the absorbed part.

    // roughness control how spread is the reflection
    // a higher roughness will cause reflect more widely
    float roughness = hit.material.rough;
    float max_angle = roughness * M_PI; // max spread out angel in radians

    srand(max_angle * reflection_intensity);

    Ray perfect_output = Reflect(input, hit.N, hit.interacts);
    float limit_angle = dot(hit.N, perfe ct_output.dir);

    int rayNum = 5; // could adjust here based on computation resource
    while (rayNum > 0) // random choose 5 ray in the range of maxangle (determined by roughness)
    {
        float randx = (float)rand() / (RAND_MAX);
        float randy = (float)rand() / (RAND_MAX);
        float randz = (float)rand() / (RAND_MAX);
        Dir3D reflect_ray(randx, randy, randz);
        reflect_ray = reflect_ray.normalized();

        float angle_to_reflect_ray = dot(reflect_ray, perfect_output.dir);
        if (angle_to_reflect_ray < max_angle) // a suitable angel
        {
            Ray output_ray(hit.interacts + hit.N * 0.01, reflect_ray);
            output_ray.intensity = 1.0 / rayNum; // all these ray share the same intensity
            result.push_back(output_ray);
            rayNum--;
        }
    }

    return result;
}



