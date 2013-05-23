#include <GL/glfw3.h>
#include <glhck/glhck.h>
#include <stdlib.h>

/***
 * Kazmath extension
 * When stuff is tested and working
 * Try to get much as we can to upstream
 ***/

typedef struct kmAABBExtent {
   kmVec3 point;
   kmVec3 extent;
} kmAABBExtent;

typedef struct kmSphere {
   kmVec3 point;
   kmScalar radius;
} kmSphere;

typedef struct kmEllipse {
   kmVec3 point;
   kmVec3 radius;
} kmEllipse;

typedef struct kmOBB {
   struct kmAABBExtent aabb;
   kmVec3 orientation[3];
} kmOBB;

typedef struct kmTriangle {
   kmVec3 v1, v2, v3;
} kmTriangle;

kmVec3* kmVec3Divide(kmVec3 *pOut, const kmVec3 *pV1, const kmVec3 *pV2)
{
   pOut->x = pV1->x / pV2->x;
   pOut->y = pV1->y / pV2->y;
   pOut->z = pV1->z / pV2->z;
   return pOut;
}

kmVec3* kmVec3Multiply(kmVec3 *pOut, const kmVec3 *pV1, const kmVec3 *pV2)
{
   pOut->x = pV1->x * pV2->x;
   pOut->y = pV1->y * pV2->y;
   pOut->z = pV1->z * pV2->z;
   return pOut;
}

static kmVec3* kmVec3Min(kmVec3 *pOut, const kmVec3 *pIn, const kmVec3 *pV1)
{
   pOut->x = min(pIn->x, pV1->x);
   pOut->y = min(pIn->y, pV1->y);
   pOut->z = min(pIn->z, pV1->z);
}

static kmVec3* kmVec3Max(kmVec3 *pOut, const kmVec3 *pIn, const kmVec3 *pV1)
{
   pOut->x = max(pIn->x, pV1->x);
   pOut->y = max(pIn->y, pV1->y);
   pOut->z = max(pIn->z, pV1->z);
}

static void kmSwap(kmScalar *a, kmScalar *b)
{
   const kmScalar tmp = *a;
   *a = *b;
   *b = tmp;
}

kmScalar kmPlaneDistanceTo(const kmPlane *pIn, const kmVec3 *pV1)
{
   const kmVec3 normal = {pIn->a, pIn->b, pIn->c};
   return kmVec3Dot(pV1, &normal) * pIn->d;
}

kmSphere* kmSphereFromAABB(kmSphere *sphere, const kmAABB *aabb)
{
   kmAABBCentre(aabb, &sphere->point);
   sphere->radius = kmAABBDiameterX(aabb);
   sphere->radius = max(sphere->radius, kmAABBDiameterY(aabb));
   sphere->radius = max(sphere->radius, kmAABBDiameterZ(aabb));
   return sphere;
}

kmBool kmSphereIntersectsSphere(const kmSphere *a, const kmSphere *b)
{
   kmVec3 vector;
   kmScalar distance, radiusSum;
   kmVec3Subtract(&vector, &a->point, &b->point);
   distance = kmVec3LengthSq(&vector);
   radiusSum = a->radius + b->radius;
   return (distance <= radiusSum * radiusSum);
}

kmBool kmOBBTestOBB(const kmOBB *a, const kmOBB *b)
{
   int i, j;
   kmMat3 mat;
   kmVec3 tmp, translation;

   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
         mat.mat[i*j] = kmVec3Dot(&a->orientation[i], &b->orientation[j]);

   kmVec3Subtract(&tmp, &a->aabb.point, &b->aabb.point);
   translation.x = kmVec3Dot(&tmp, &a->orientation[0]);
   translation.y = kmVec3Dot(&tmp, &a->orientation[2]);
   translation.z = kmVec3Dot(&tmp, &a->orientation[2]);
}

kmBool kmAABBIntersectsAABB(const kmAABB *a, const kmAABB *b)
{
   if (a->max.x < b->min.x || a->min.x > b->max.x) return KM_FALSE;
   if (a->max.y < b->min.y || a->min.y > b->max.y) return KM_FALSE;
   if (a->max.z < b->min.z || a->min.z > b->max.z) return KM_FALSE;
   return KM_TRUE;
}

kmBool kmAABBExtentIntersectsAABBExtent(const kmAABBExtent *a, const kmAABBExtent *b)
{
   if (abs(a->point.x - b->point.x) > (a->extent.x + b->extent.x)) return KM_FALSE;
   if (abs(a->point.y - b->point.y) > (a->extent.y + b->extent.y)) return KM_FALSE;
   if (abs(a->point.z - b->point.z) > (a->extent.z + b->extent.z)) return KM_FALSE;
   return KM_TRUE;
}

/***
 * GLHCK kazmath integration
 * Create glhck objects from kazmath primitives
 * Move to glhck when tested and working
 ***/

glhckObject* glhckCubeFromKazmathAABB(const kmAABB *aabb)
{
   kmVec3 center;
   glhckObject *o = glhckCubeNew(1.0);
   kmAABBCentre(aabb, &center);
   glhckObjectPosition(o, &center);
   glhckObjectScalef(o,
         kmAABBDiameterX(aabb)*0.5,
         kmAABBDiameterY(aabb)*0.5,
         kmAABBDiameterZ(aabb)*0.5);
   return o;
}

glhckObject* glhckCubeFromKazmathAABBExtent(const kmAABBExtent *aabb)
{
   glhckObject *o = glhckCubeNew(1.0);
   glhckObjectPosition(o, &aabb->point);
   glhckObjectScalef(o, aabb->extent.x, aabb->extent.y, aabb->extent.z);
   return o;
}

static const int WIDTH = 800;
static const int HEIGHT = 480;
static char RUNNING = 1;
static void windowCloseCallback(GLFWwindow *window)
{
   (void)window;
   RUNNING = 0;
}

typedef struct AABBvsAABB {
   const kmAABB *a, *b;
   glhckObject *oa, *ob;
} AABBvsAABB;

typedef struct AABBEvsAABBE {
   const kmAABBExtent *a, *b;
} AABBEvsAABBE;

typedef struct OBBvsOBB {
   const kmOBB *a, *b;
} OBBvsOBB;

typedef struct SphereVsSphere {
   const kmSphere *a, *b;
} SphereVsSphere;

static void textToObject(glhckText *text, unsigned int font, int size, glhckObject *object, const char *str)
{
   const kmVec3 *posit = glhckObjectGetPosition(object);
   const kmVec3 *scale = glhckObjectGetScale(object);
   glhckTextStash(text, font, size, posit->x-scale->x, posit->y+scale->y+size, str, NULL);
   glhckTextRender(text);
   glhckTextClear(text);
}

static void run(GLFWwindow *window)
{
   unsigned int i, textSize, font;
   char keyLock = 0;
   glhckText *text;
   kmMat4 ortho;

   text = glhckTextNew(1024, 1024);
   font = glhckTextFontNewKakwafont(text, &textSize);

   #define LENGTH(X) (sizeof X / sizeof X[0])
   AABBvsAABB aabbTests[] = {
      {.a = (&(kmAABB){.min = {WIDTH*0.2-50,HEIGHT/2-50,0}, .max = {WIDTH*0.2+50,HEIGHT/2+50,0}}),
       .b = (&(kmAABB){.min = {WIDTH*0.8-50,HEIGHT/2-50,0}, .max = {WIDTH*0.8+50,HEIGHT/2+50,0}})},
      {.a = (&(kmAABB){.min = {WIDTH*0.2-50,HEIGHT/2-50,0}, .max = {WIDTH*0.2+50,HEIGHT/2+50,0}}),
       .b = (&(kmAABB){.min = {WIDTH*0.4-50,HEIGHT/2-50,0}, .max = {WIDTH*0.8+50,HEIGHT/2+50,0}})},
      {.a = (&(kmAABB){.min = {WIDTH*0.2-50,HEIGHT/2-50,0}, .max = {WIDTH*0.2+50,HEIGHT/2+50,0}}),
       .b = (&(kmAABB){.min = {WIDTH*0.2+00,HEIGHT/2+20,0}, .max = {WIDTH*0.2+100,HEIGHT/2+120,0}})},
   };
   SphereVsSphere sphereTests[] = {
      {.a = (&(kmSphere){.point = {0,0,0}, .radius = 5}),
       .b = (&(kmSphere){.point = 6, .radius = 2})},
   };

   for (i = 0; i < LENGTH(aabbTests); ++i) {
      aabbTests[i].oa = glhckCubeFromKazmathAABB(aabbTests[i].a);
      aabbTests[i].ob = glhckCubeFromKazmathAABB(aabbTests[i].b);
      glhckObjectDrawAABB(aabbTests[i].oa, 1);
      glhckObjectDrawAABB(aabbTests[i].ob, 1);
      glhckMaterial *amat = glhckMaterialNew(NULL);
      glhckMaterial *bmat = glhckMaterialNew(NULL);
      glhckObjectMaterial(aabbTests[i].oa, amat);
      glhckObjectMaterial(aabbTests[i].ob, bmat);
   }

   /* lets use 2D projection for easier visualization of data
    * the text is left-handed, the world in glhck is left-handed, but Y is reversed (positive == UP)
    * to get same ortho for world as text, we do this
    *
    *  TEXT    WORLD/3D
    *   -y        y
    * -x | x   -x | x
    *    y       -y
    *
    * */
   kmMat4OrthographicProjection(&ortho, 0, WIDTH, HEIGHT, 0, -1, 1);
   glhckRenderProjectionOnly(&ortho);

   i = 0;
   while(RUNNING && glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS) {
      glfwPollEvents();

      if (glfwGetKey(window, GLFW_KEY_SPACE)) {
         if (!keyLock && ++i >= LENGTH(aabbTests)) i = 0;
         keyLock = 1;
      } else {
         keyLock = 0;
      }

      if (i < LENGTH(aabbTests)) {
         glhckMaterial *amat = glhckObjectGetMaterial(aabbTests[i].oa);
         glhckMaterial *bmat = glhckObjectGetMaterial(aabbTests[i].ob);

         if (!kmAABBIntersectsAABB(aabbTests[i].a, aabbTests[i].b)) {
            textToObject(text, font, textSize, aabbTests[i].oa, "[A] No intersection");
            textToObject(text, font, textSize, aabbTests[i].ob, "[B] No intersection");
         } else {
            textToObject(text, font, textSize, aabbTests[i].oa, "[A] Intersection");
            textToObject(text, font, textSize, aabbTests[i].ob, "[B] Intersection");
            glhckMaterialDiffuseb(amat, 255, 0, 0, 255);
            glhckMaterialDiffuseb(bmat, 255, 0, 0, 255);
         }

         glhckObjectRender(aabbTests[i].oa);
         glhckObjectRender(aabbTests[i].ob);
         glhckMaterialDiffuseb(amat, 0, 0, 255, 255);
         glhckMaterialDiffuseb(bmat, 0, 255, 0, 255);
      }

      glhckTextStash(text, font, textSize, 0, textSize, "colhck - coltest.c", NULL);
      glhckTextStash(text, font, textSize, 0, textSize*2, "press [space] for next test case", NULL);
      glhckTextRender(text);
      glhckTextClear(text);

      glfwSwapBuffers(window);
      glhckRenderClear(GLHCK_DEPTH_BUFFER | GLHCK_COLOR_BUFFER);
   }
}

int main(int argc, char **argv)
{
   GLFWwindow *window;
   if (!glfwInit())
      return EXIT_FAILURE;

   glfwWindowHint(GLFW_DEPTH_BITS, 24);
   if (!(window = glfwCreateWindow(WIDTH, HEIGHT, "colhck - coltest.c", NULL, NULL)))
      return EXIT_FAILURE;

   glfwMakeContextCurrent(window);
   glfwSwapInterval(0);

   if (!glhckContextCreate(argc, argv))
      return EXIT_FAILURE;

   if (!glhckDisplayCreate(WIDTH, HEIGHT, GLHCK_RENDER_OPENGL))
      return EXIT_FAILURE;

   glfwSetWindowCloseCallback(window, windowCloseCallback);
   run(window);

   glhckContextTerminate();
   glfwTerminate();
   return EXIT_SUCCESS;
}
