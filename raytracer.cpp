#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <climits>

using namespace std;
using namespace parser;

typedef unsigned char RGB[3];

class Ray
{
public:
    Vec3f origin{};
    Vec3f direction{};
    bool isShadow;

    Ray(Vec3f origin, Vec3f direction, bool isShadow = false){
        this->origin = origin;
        this->direction = direction;
        this-> isShadow = isShadow;
    }
};

class Hit
{
public:
    bool hit; //returns true if ray hits an object.
    float t;
    Vec3f x; //Intersection point.
    Vec3f normal; // surface normal (normalized)
    int materialID;
};


Hit SphereIntersection(Ray ray, Sphere sphere);

Vec3f CrossProduct(Vec3f first, Vec3f second)
{
    Vec3f result = {0,0,0};
    result.x = first.y * second.z - first.z * second.y;
    result.y = first.z * second.x - first.x * second.z;
    result.z = first.x * second.y - first.y * second.x;

    return result;
}



float DotProduct(Vec3f first, Vec3f second){
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

float Determinant(Vec3f v1, Vec3f v2, Vec3f v3)
{
    return v1.x * (v2.y * v3.z - v3.y * v2.z) - v1.y * (v2.x * v3.z - v3.x * v2.z) + v1.z * (v2.x * v3.y - v3.x * v2.y);
}

Vec3f NegateVector(Vec3f vector)
{
    return {-vector.x, -vector.y, -vector.z};
}

Vec3f SumVectors(Vec3f first, Vec3f second)
{
    return {first.x + second.x, first.y + second.y, first.z + second.z };
}

Vec3f Normalize(Vec3f vector)
{
    float magnitude;
    magnitude = sqrtf(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2));
    return {vector.x / magnitude, vector.y / magnitude, vector.z / magnitude };
}

Vec3f SubtractVectors(Vec3f first, Vec3f second)
{
    Vec3f result = {first.x - second.x, first.y - second.y, first.z - second.z };
    return Normalize(result);

}


Vec3f FindNormal(Vec3f a, Vec3f b, Vec3f c)
{
    return Normalize(CrossProduct(SubtractVectors(b,a), SubtractVectors(c,a)));  // May not be true. The direction may be different. need to check.
}

float PixelHeight(const Camera& camera){
    // x-> left y-> right z-> bottom w-> top
    return (camera.near_plane.w - camera.near_plane.z)/(float) camera.image_height;
}

float PixelWidth(const Camera& camera){
    return (camera.near_plane.y - camera.near_plane.x)/(float) camera.image_width;
}

Vec3f MultiplyVector(Vec3f v, float c)
{
    return {v.x * c, v.y * c, v.z * c};
}


Vec3f** BuildPixelMatrix(Camera camera)
{
    float pw = PixelWidth(camera);
    float ph = PixelHeight(camera);

    Vec3f u = CrossProduct(camera.up, NegateVector(camera.gaze));
    // -----Initialize the matrix ----------
    Vec3f** result = new Vec3f*[camera.image_width];
    for (int i =0; i< camera.image_width; i++)
    {
        result[i] = new Vec3f[camera.image_height];
    }


    for(int i=0; i < camera.image_width; i++)
    {
        for(int j=0; j < camera.image_height; j++)
        {
            result[i][j].x = camera.position.x +
                                u.x * (camera.near_plane.x + ( i + 0.5f) * pw )+
                                camera.up.x * (camera.near_plane.w - (j + 0.5f)*ph)+
                                camera.gaze.x * camera.near_distance;
            result[i][j].y = camera.position.y +
                                u.y * (camera.near_plane.x + ( i + 0.5f) * pw )+
                                camera.up.y * (camera.near_plane.w - (j + 0.5f)*ph)+
                                camera.gaze.y * camera.near_distance;
            result[i][j].z = camera.position.z +
                                u.z * (camera.near_plane.x + ( i + 0.5f) * pw )+
                                camera.up.z * (camera.near_plane.w - (j + 0.5f) * ph)+
                                camera.gaze.z * camera.near_distance;
        }



    }
    return result;
}

Hit TriangleIntersection(Ray ray, Vec3f a, Vec3f b, Vec3f c, int materialID, Vec3f normal)
{
    Hit result_hit{};
    float alpha, beta, gamma, t, detA;
    Vec3f minusD = NegateVector(ray.direction);
    Vec3f a_minus_c = SubtractVectors(a,c);
    Vec3f b_minus_c = SubtractVectors(b,c);
    Vec3f o_minus_c = SubtractVectors(ray.origin,c);
    Vec3f alphaMatrix[3];
    Vec3f betaMatrix[3];
    Vec3f tMatrix[3];
    // ---Calculate Determinant of A---
    detA = Determinant(a_minus_c, b_minus_c, minusD);
    if(detA == 0.0)
    {
        result_hit.hit = false;
        return result_hit;
    }

    alpha = Determinant(o_minus_c, b_minus_c, minusD) / detA;
    beta = Determinant(a_minus_c, o_minus_c, minusD) / detA;
    t = Determinant(a_minus_c, b_minus_c, o_minus_c) / detA;
    gamma = 1.0 - (alpha + beta);

    if((0 <= alpha && alpha <= 1) && (0 <= beta && beta <= 1) && (0 <= gamma && gamma <= 1))
    {
        Vec3f intersectionPoint = SumVectors(ray.origin, MultiplyVector(ray.direction, t));
        result_hit.hit = true;
        result_hit.t = t;
        result_hit.x = intersectionPoint;
        result_hit.normal = normal;
        result_hit.materialID = materialID;
        return result_hit;
    }

    result_hit.hit = false;
    return result_hit;
}

Hit SphereIntersection(Ray ray, Sphere sphere, Vec3f center)
{
    Hit result_hit{};
    float A,B,C ,t;
    Vec3f o_minus_c = SubtractVectors(ray.origin, center);
    A = DotProduct(ray.direction,ray.direction);
    B = 2 * DotProduct(ray.direction, o_minus_c);
    C = DotProduct(o_minus_c, o_minus_c) - sphere.radius * sphere.radius;

    float delta = B*B - 4 * A * C ;
    if(delta < 0.0)
    {
        result_hit.hit = false;
        return result_hit;
    }
    if(delta == 0.0)
    {
        t = (-B)/(2 * A);
        Vec3f intersectionPoint = SumVectors(ray.origin, MultiplyVector(ray.direction, t));
        result_hit.hit = true,
        result_hit.t = t;
        result_hit.materialID = sphere.material_id;
        result_hit.x = intersectionPoint;
        result_hit.normal = SubtractVectors(intersectionPoint, center);
        return result_hit;
    }
    float t1 = (-B + sqrtf(delta))/(2 * A);
    float t2 = (-B - sqrtf(delta))/(2 * A);
    if(t1 < 0 && t2 < 0)
    {
        result_hit.hit = false;
        return result_hit;
    }
    if(t1 > 0 && t2 > 0)    t = fmin(t1,t2);
    else    t = fmax(t1,t2);

    Vec3f intersectionPoint = SumVectors(ray.origin, MultiplyVector(ray.direction, t));
    result_hit.hit = true,
    result_hit.t = t;
    result_hit.materialID = sphere.material_id;
    result_hit.x = intersectionPoint;
    result_hit.normal = SubtractVectors(intersectionPoint, center);

    return result_hit;
}








int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;
    scene.loadFromXml(argv[1]);
    int numCameras = scene.cameras.size();

    //--- Calculate the surface normal of each triangle and sphere beforehand.
    int triangleSize = scene.triangles.size();
    int sphereSize = scene.spheres.size();
    Vec3f triangleNormals[triangleSize];
    for(int i = 0; i < triangleSize; i++)
    {
        Vec3f a = scene.vertex_data[scene.triangles[i].indices.v0_id - 1];
        Vec3f b = scene.vertex_data[scene.triangles[i].indices.v1_id - 1];
        Vec3f c = scene.vertex_data[scene.triangles[i].indices.v2_id - 1];

        triangleNormals[i] = FindNormal(a,b,c);
    }




    for(int cameraIndex = 0; cameraIndex < numCameras; cameraIndex++)
    {

        Camera currentCamera = scene.cameras[i];
        int image_width = currentCamera.image_width;
        int image_height = currentCamera.image_height;
        unsigned char* image = new unsigned char [image_width * image_height * 3];
        Vec3f** pixelMatrix = BuildPixelMatrix(currentCamera);

        for(int i = 0; i < image_width; i++ )
        {
            for(int j= 0; j < image_height; j++)
            {
                //Compute the ray
                Ray ray = Ray(currentCamera.position, Normalize(SubtractVectors(pixelMatrix[i][j], currentCamera.position)));

                float t_min = numeric_limits<float>::infinity();
                Hit hit{};
                // ---Intersect with Triangles.
                for(int t_index = 0; t_index < triangleSize; t_index++)
                {

                    Triangle triangle = scene.triangles[t_index];
                    Vec3f a = scene.vertex_data[triangle.indices.v0_id - 1];
                    Vec3f b = scene.vertex_data[triangle.indices.v1_id - 1];
                    Vec3f c = scene.vertex_data[triangle.indices.v2_id - 1];
                    float t;

                    Hit new_hit = TriangleIntersection(ray, a, b, c, triangle.material_id, triangleNormals[t_index]);
                    t = new_hit.t;
                    if(new_hit.hit)
                    {
                        if (t < t_min)
                        {
                            t_min = t;
                            hit = new_hit;
                        }
                    }
                }

                // ---Intersect with Spheres
                for(int s_index = 0; s_index < sphereSize; s_index++)
                {
                    Sphere sphere = scene.spheres[s_index];
                    float t;
                    Vec3f center = scene.vertex_data[sphere.center_vertex_id - 1];
                    Hit new_hit = SphereIntersection(ray, sphere, center);
                    t = new_hit.t;

                    if(new_hit.hit)
                    {
                        if (t < t_min)
                        {
                            t_min = t;
                            hit = new_hit;
                        }
                    }
                }
                // TODO






            }
        }

    }
    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };

//    int width = 640, height = 480;
//    int columnWidth = width / 8;



    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }
    if (i == 1)
    write_ppm("test.ppm", image, width, height);

}


