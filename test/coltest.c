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

kmMat3* kmOBBGetMat3(const kmOBB *pIn, kmMat3 *pOut)
{
   int i;
   for (i = 0; i < 3; ++i) {
      pOut->mat[i*3+0] = pIn->orientation[i].x;
      pOut->mat[i*3+1] = pIn->orientation[i].y;
      pOut->mat[i*3+2] = pIn->orientation[i].z;
   }
   return pOut;
}

kmBool kmOBBIntersectsOBB(const kmOBB *a, const kmOBB *b)
{
   int i, j;
   kmScalar ra, rb;
   kmMat3 mat, absMat;
   kmVec3 tmp, translation;

   /* compute rotation matrix expressing b in a's coordinate frame */
   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
         mat.mat[i*3+j] = kmVec3Dot(&a->orientation[i], &b->orientation[j]);

   /* bring translations into a's coordinate frame */
   kmVec3Subtract(&tmp, &a->aabb.point, &b->aabb.point);
   translation.x = kmVec3Dot(&tmp, &a->orientation[0]);
   translation.y = kmVec3Dot(&tmp, &a->orientation[2]);
   translation.z = kmVec3Dot(&tmp, &a->orientation[2]);

   /* compute common subexpressions. add in and epsilon term to
    * counteract arithmetic errors when two edges are parallel and
    * their cross product is (near) null. */
   for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
         absMat.mat[i*j] = abs(mat.mat[i*j]) + kmEpsilon;

   /* test axes L = A0, L = A1, L = A2 */
   for (i = 0; i < 3; ++i) {
      ra = (i==0?a->aabb.extent.x:i==1?a->aabb.extent.y:a->aabb.extent.z);
      rb = b->aabb.extent.x * absMat.mat[i*3+0] + b->aabb.extent.y * absMat.mat[i*3+1] + b->aabb.extent.z * absMat.mat[i*3+2];
      if (abs((i==0?translation.x:i==1?translation.y:translation.z)) > ra + rb) return KM_FALSE;
   }

   /* test axes L = B0, L = B1, L = B2 */
   for (i = 0; i < 3; ++i) {
      ra = a->aabb.extent.x * absMat.mat[0*3+i] + a->aabb.extent.y * absMat.mat[1*3+i] + a->aabb.extent.z * absMat.mat[2*3+i];
      rb = (i==0?b->aabb.extent.x:i==1?b->aabb.extent.y:b->aabb.extent.z);
      if (abs(translation.x * mat.mat[0*3+i] + translation.y * mat.mat[1*3+i] + translation.z * mat.mat[2*3+i]) > ra + rb) return KM_FALSE;
   }

   /* test axis L = A0 x B0 */
   ra = a->aabb.extent.y * absMat.mat[2*3+0] + a->aabb.extent.z * absMat.mat[1*3+0];
   rb = b->aabb.extent.y * absMat.mat[0*3+2] + b->aabb.extent.z * absMat.mat[0*3+1];
   if (abs(translation.z * mat.mat[1*3+0] - translation.y * mat.mat[2*3+0]) > ra + rb) return KM_FALSE;

   /* test axis L = A0 x B1 */
   ra = a->aabb.extent.y * absMat.mat[2*3+1] + a->aabb.extent.z * absMat.mat[1*3+1];
   rb = b->aabb.extent.x * absMat.mat[0*3+2] + b->aabb.extent.z * absMat.mat[0*3+0];
   if (abs(translation.z * mat.mat[1*3+1] - translation.y * mat.mat[2*3+1]) > ra + rb) return KM_FALSE;

   /* test axis L = A0 x B2 */
   ra = a->aabb.extent.y * absMat.mat[2*3+2] + a->aabb.extent.z * absMat.mat[1*3+2];
   rb = b->aabb.extent.x * absMat.mat[0*3+1] + b->aabb.extent.y * absMat.mat[0*3+0];
   if (abs(translation.z * mat.mat[1*3+2] - translation.y * mat.mat[2*3+2]) > ra + rb) return KM_FALSE;

   /* test axis L = A1 x B0 */
   ra = a->aabb.extent.x * absMat.mat[2*3+0] + a->aabb.extent.z * absMat.mat[0*3+0];
   rb = b->aabb.extent.y * absMat.mat[1*3+2] + b->aabb.extent.z * absMat.mat[1*3+1];
   if (abs(translation.x * mat.mat[2*3+0] - translation.z * mat.mat[0*3+0]) > ra + rb) return KM_FALSE;

   /* test axis L = A1 x B1 */
   ra = a->aabb.extent.x * absMat.mat[2*3+1] + a->aabb.extent.z * absMat.mat[0*3+1];
   rb = b->aabb.extent.x * absMat.mat[1*3+2] + b->aabb.extent.z * absMat.mat[1*3+0];
   if (abs(translation.x * mat.mat[2*3+1] - translation.z * mat.mat[0*3+1]) > ra + rb) return KM_FALSE;

   /* test axis L = A1 x B2 */
   ra = a->aabb.extent.x * absMat.mat[2*3+2] + a->aabb.extent.z * absMat.mat[0*3+2];
   rb = b->aabb.extent.x * absMat.mat[1*3+1] + b->aabb.extent.y * absMat.mat[1*3+0];
   if (abs(translation.x * mat.mat[2*3+2] - translation.z * mat.mat[0*3+2]) > ra + rb) return KM_FALSE;

   /* test axis L = A2 x B0 */
   ra = a->aabb.extent.x * absMat.mat[1*3+0] + a->aabb.extent.y * absMat.mat[0*3+0];
   rb = b->aabb.extent.y * absMat.mat[2*3+2] + b->aabb.extent.z * absMat.mat[2*3+1];
   if (abs(translation.y * mat.mat[0*3+0] - translation.x * mat.mat[1*3+0]) > ra + rb) return KM_FALSE;

   /* test axis L = A2 x B1 */
   ra = a->aabb.extent.x * absMat.mat[1*3+1] + a->aabb.extent.y * absMat.mat[0*3+1];
   rb = b->aabb.extent.x * absMat.mat[2*3+2] + b->aabb.extent.z * absMat.mat[2*3+0];
   if (abs(translation.y * mat.mat[0*3+1] - translation.x * mat.mat[1*3+1]) > ra + rb) return KM_FALSE;

   /* test axis L = A2 x B2 */
   ra = a->aabb.extent.x * absMat.mat[1*3+2] + a->aabb.extent.y * absMat.mat[0*3+2];
   rb = b->aabb.extent.x * absMat.mat[2*3+1] + b->aabb.extent.y * absMat.mat[2*3+0];
   if (abs(translation.y * mat.mat[0*3+2] - translation.x * mat.mat[1*3+2]) > ra + rb) return KM_FALSE;

   /* no seperating axis found */
   return KM_TRUE;
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

glhckObject* glhckCubeFromKazmathOBB(const kmOBB *obb)
{
   kmMat3 mat;
   kmVec3 rot = {1,1,1};
   glhckObject *o = glhckCubeFromKazmathAABBExtent(&obb->aabb);
   kmOBBGetMat3(obb, &mat);
   kmVec3MultiplyMat3(&rot, &rot, &mat);
   glhckObjectRotation(o, &rot);
   return o;
}

glhckObject* glhckCubeFromKazmathSphere(const kmSphere *sphere)
{
   glhckObject *o = glhckSphereNew(sphere->radius);
   glhckObjectPosition(o, &sphere->point);
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
   char intersection;
} AABBvsAABB;

typedef struct AABBEvsAABBE {
   const kmAABBExtent *a, *b;
   glhckObject *oa, *ob;
   char intersection;
} AABBEvsAABBE;

typedef struct OBBvsOBB {
   const kmOBB *a, *b;
   glhckObject *oa, *ob;
   char intersection;
} OBBvsOBB;

typedef struct SphereVsSphere {
   const kmSphere *a, *b;
   glhckObject *oa, *ob;
   char intersection;
} SphereVsSphere;

static void textToObject(glhckText *text, unsigned int font, int size, glhckObject *object, const char *str)
{
   const kmVec3 *posit = glhckObjectGetPosition(object);
   const kmAABB *aabb = glhckObjectGetAABB(object);
   glhckTextStash(text, font, size, aabb->min.x, aabb->max.y+1+size, str, NULL);
   glhckTextRender(text);
   glhckTextClear(text);
}

static void run(GLFWwindow *window)
{
   typedef enum primitiveType {
      PRIMITIVE_AABB,
      PRIMITIVE_AABBE,
      PRIMITIVE_OBB,
      PRIMITIVE_SPHERE,
      PRIMITIVE_LAST,
   } primitiveType;

   primitiveType primitive = PRIMITIVE_AABB;
   unsigned int i, textSize, font, testsInPrimitive = 0;
   char keyLock = 0, expectIntersection = 0, didIntersect = 0;
   char *primitiveStr = NULL;
   glhckText *text;
   kmMat4 ortho;

   text = glhckTextNew(1024, 1024);
   font = glhckTextFontNewKakwafont(text, &textSize);

   #define LENGTH(X) (sizeof X / sizeof X[0])
   AABBvsAABB aabbTests[] = {
      {.a = (&(kmAABB){.min = {WIDTH*0.2-50,HEIGHT/2-50,0}, .max = {WIDTH*0.2+50,HEIGHT/2+50,0}}),
       .b = (&(kmAABB){.min = {WIDTH*0.8-50,HEIGHT/2-50,0}, .max = {WIDTH*0.8+50,HEIGHT/2+50,0}}),
       .intersection = 0},
      {.a = (&(kmAABB){.min = {WIDTH*0.2-50,HEIGHT/2-50,0}, .max = {WIDTH*0.2+50,HEIGHT/2+50,0}}),
       .b = (&(kmAABB){.min = {WIDTH*0.4-50,HEIGHT/2-50,0}, .max = {WIDTH*0.8+50,HEIGHT/2+50,0}}),
       .intersection = 0},
      {.a = (&(kmAABB){.min = {WIDTH*0.2-50,HEIGHT/2-50,0}, .max = {WIDTH*0.2+50,HEIGHT/2+50,0}}),
       .b = (&(kmAABB){.min = {WIDTH*0.2+00,HEIGHT/2+20,0}, .max = {WIDTH*0.2+100,HEIGHT/2+120,0}}),
       .intersection = 1},
   };
   AABBEvsAABBE aabbeTests[] = {
      {.a = (&(kmAABBExtent){.point = {WIDTH*0.2,HEIGHT/2,0}, .extent = {50,50,0}}),
       .b = (&(kmAABBExtent){.point = {WIDTH*0.8,HEIGHT/2,0}, .extent = {50,50,0}}),
       .intersection = 0},
      {.a = (&(kmAABBExtent){.point = {WIDTH*0.2,HEIGHT/2,0}, .extent = {50,50,0}}),
       .b = (&(kmAABBExtent){.point = {WIDTH*0.6,HEIGHT/2,0}, .extent = {200,50,0}}),
       .intersection = 0},
      {.a = (&(kmAABBExtent){.point = {WIDTH*0.2,HEIGHT/2,0}, .extent = {50,50,0}}),
       .b = (&(kmAABBExtent){.point = {WIDTH*0.2+50,HEIGHT/2+50,0}, .extent = {50,50,0}}),
       .intersection = 1},
   };
   OBBvsOBB obbTests[] = {
      {.a = (&(kmOBB){.aabb = *aabbeTests[0].a,
            .orientation[0] = {1,1,45},
            .orientation[1] = {1,1,45},
            .orientation[2] = {1,1,45}}),
       .b = (&(kmOBB){.aabb = *aabbeTests[0].b,
            .orientation[0] = {1,1,45},
            .orientation[1] = {1,1,45},
            .orientation[2] = {1,1,45}}),
      .intersection = 0},
      {.a = (&(kmOBB){.aabb = *aabbeTests[1].a,
            .orientation[0] = {1,1,45},
            .orientation[1] = {1,1,45},
            .orientation[2] = {1,1,45}}),
       .b = (&(kmOBB){.aabb = *aabbeTests[1].b,
            .orientation[0] = {1,1,45},
            .orientation[1] = {1,1,45},
            .orientation[2] = {1,1,45}}),
      .intersection = 0},
      {.a = (&(kmOBB){.aabb = *aabbeTests[2].a,
            .orientation[0] = {1,1,45},
            .orientation[1] = {1,1,45},
            .orientation[2] = {1,1,45}}),
       .b = (&(kmOBB){.aabb = *aabbeTests[2].b,
            .orientation[0] = {1,1,45},
            .orientation[1] = {1,1,45},
            .orientation[2] = {1,1,45}}),
      .intersection = 1},
   };
   SphereVsSphere sphereTests[] = {
      {.a = (&(kmSphere){.point = {WIDTH*0.2,HEIGHT/2,0}, .radius = 50}),
       .b = (&(kmSphere){.point = {WIDTH*0.8,HEIGHT/2,0}, .radius = 50}),
       .intersection = 0},
      {.a = (&(kmSphere){.point = {WIDTH*0.2,HEIGHT/2,0}, .radius = 50}),
       .b = (&(kmSphere){.point = {WIDTH*0.42,HEIGHT/2,0}, .radius = 120}),
       .intersection = 0},
      {.a = (&(kmSphere){.point = {WIDTH*0.2,HEIGHT/2,0}, .radius = 50}),
       .b = (&(kmSphere){.point = {WIDTH*0.2+50,HEIGHT/2+50,0}, .radius = 50}),
       .intersection = 1},
   };

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
   kmMat4OrthographicProjection(&ortho, 0, WIDTH, HEIGHT, 0, -500, 500);
   glhckRenderProjectionOnly(&ortho);

   #define createPrimitives(i, primitives, creatFunc) \
      for (i = 0; i < LENGTH(primitives); ++i) { \
         primitives[i].oa = creatFunc(primitives[i].a); \
         primitives[i].ob = creatFunc(primitives[i].b); \
         glhckObjectDrawAABB(primitives[i].oa, 1); \
         glhckObjectDrawAABB(primitives[i].ob, 1); \
         glhckMaterial *amat = glhckMaterialNew(NULL); \
         glhckMaterial *bmat = glhckMaterialNew(NULL); \
         glhckObjectMaterial(primitives[i].oa, amat); \
         glhckObjectMaterial(primitives[i].ob, bmat); \
      }

   #define testPrimitive(i, primitives, interFunc, testsInPrimitive, expectIntersection, didIntersect) \
      if (i < (testsInPrimitive = LENGTH(primitives))) { \
         glhckMaterial *amat = glhckObjectGetMaterial(primitives[i].oa); \
         glhckMaterial *bmat = glhckObjectGetMaterial(primitives[i].ob); \
         if (!interFunc(primitives[i].a, primitives[i].b)) { \
            textToObject(text, font, textSize, primitives[i].oa, "[A] No intersection"); \
            textToObject(text, font, textSize, primitives[i].ob, "[B] No intersection"); \
            didIntersect = 0; \
         } else { \
            textToObject(text, font, textSize, primitives[i].oa, "[A] Intersection"); \
            textToObject(text, font, textSize, primitives[i].ob, "[B] Intersection"); \
            glhckMaterialDiffuseb(amat, 255, 0, 0, 255); \
            glhckMaterialDiffuseb(bmat, 255, 0, 0, 255); \
            didIntersect = 1; \
         } \
         glhckObjectRender(primitives[i].oa); \
         glhckObjectRender(primitives[i].ob); \
         glhckMaterialDiffuseb(amat, 0, 0, 255, 255); \
         glhckMaterialDiffuseb(bmat, 0, 255, 0, 255); \
         expectIntersection = primitives[i].intersection; \
      }

   /* create test primitives */
   createPrimitives(i, aabbTests, glhckCubeFromKazmathAABB);
   createPrimitives(i, aabbeTests, glhckCubeFromKazmathAABBExtent);
   createPrimitives(i, obbTests, glhckCubeFromKazmathOBB);
   createPrimitives(i, sphereTests, glhckCubeFromKazmathSphere);

   i = 0;
   while(RUNNING && glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS) {
      glfwPollEvents();

      /* test intersection */
      switch (primitive) {
         case PRIMITIVE_AABB:
            testPrimitive(i, aabbTests, kmAABBIntersectsAABB, testsInPrimitive, expectIntersection, didIntersect);
            primitiveStr = "AABBvsAABB";
            break;
         case PRIMITIVE_AABBE:
            testPrimitive(i, aabbeTests, kmAABBExtentIntersectsAABBExtent, testsInPrimitive, expectIntersection, didIntersect);
            primitiveStr = "AABBEvsAABBE";
            break;
         case PRIMITIVE_OBB:
            testPrimitive(i, obbTests, kmOBBIntersectsOBB, testsInPrimitive, expectIntersection, didIntersect);
            primitiveStr = "OBBvsOBB";
            break;
         case PRIMITIVE_SPHERE:
            testPrimitive(i, sphereTests, kmSphereIntersectsSphere, testsInPrimitive, expectIntersection, didIntersect);
            primitiveStr = "SPHEREvsSPHERE";
            break;
         default:break;
      }

      /* next intersection test */
      if (glfwGetKey(window, GLFW_KEY_SPACE)) {
         if (!keyLock && ++i >= testsInPrimitive) {
            if (++primitive >= PRIMITIVE_LAST) primitive = PRIMITIVE_AABB;
            i = 0;
         }
         keyLock = 1;
      } else {
         keyLock = 0;
      }

      /* draw hud */
      glhckTextStash(text, font, textSize, 0, textSize, "colhck - coltest.c", NULL);
      glhckTextStash(text, font, textSize, 0, textSize*2, "press [space] for next test case", NULL);
      glhckTextStash(text, font, textSize, 0, textSize*4, primitiveStr, NULL);
      if (expectIntersection) {
         glhckTextStash(text, font, textSize, 0, textSize*5, "expect intersection: YES", NULL);
      } else {
         glhckTextStash(text, font, textSize, 0, textSize*5, "expect intersection: NO", NULL);
      }
      glhckTextRender(text);
      glhckTextClear(text);

      if (didIntersect == expectIntersection) {
         glhckTextColorb(text, 0, 255, 0, 255);
         glhckTextStash(text, font, textSize, 0, textSize*6, "CORRECT", NULL);
      } else {
         glhckTextColorb(text, 255, 0, 0, 255);
         glhckTextStash(text, font, textSize, 0, textSize*6, "WRONG", NULL);
      }
      glhckTextRender(text);
      glhckTextClear(text);
      glhckTextColorb(text, 255, 255, 255, 255);

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
