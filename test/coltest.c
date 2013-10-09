#include <GLFW/glfw3.h>
#include <glhck/glhck.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define IFDO(f, x) { if (x) f(x); x = NULL; }

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

typedef struct kmCapsule {
   kmVec3 pointA;
   kmVec3 pointB;
   kmScalar radius;
} kmCapsule;

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

kmAABBExtent* kmAABBToAABBExtent(kmAABBExtent* pOut, const kmAABB* aabb)
{
  kmVec3Subtract(&pOut->extent, &aabb->max, &aabb->min);
  kmVec3Scale(&pOut->extent, &pOut->extent, 0.5f);
  kmVec3Add(&pOut->point, &aabb->min, &pOut->extent);
  return pOut;
}

kmAABB* kmAABBExtentToAABB(kmAABB* pOut, const kmAABBExtent* aabbExtent)
{
  kmVec3Subtract(&pOut->min, &aabbExtent->point, &aabbExtent->extent);
  kmVec3Add(&pOut->max, &aabbExtent->point, &aabbExtent->extent);
  return pOut;
}

static void kmMat3Print(kmMat3 *pIn)
{
   int j;
   for (j = 0; j < 3; ++j) printf("[%.2f, %.2f, %.2f]\n", pIn->mat[j*3+0], pIn->mat[j*3+1], pIn->mat[j*3+2]);
}

kmVec3* kmVec3Abs(kmVec3 *pOut, const kmVec3 *pV1)
{
   pOut->x = fabs(pV1->x);
   pOut->y = fabs(pV1->y);
   pOut->z = fabs(pV1->z);
   return pOut;
}

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

/* compute distance between segment and point */
kmScalar kmVec3LengthSqSegment(const kmVec3 *a, const kmVec3 *b, const kmVec3 *c)
{
   kmVec3 ab, ac, bc;
   kmScalar e, f;
   kmVec3Subtract(&ab, b, a);
   kmVec3Subtract(&ac, c, a);
   kmVec3Subtract(&bc, c, b);
   e = kmVec3Dot(&ac, &ab);

   /* handle cases where c projects outside ab */
   if (e <= 0.0f) return kmVec3LengthSq(&ac);

   f = kmVec3LengthSq(&ab);
   if (e >= f) return kmVec3LengthSq(&bc);

   /* handle case where c projects onto ab */
   return kmVec3LengthSq(&ac) - e * e / f;
}

kmBool kmAABBExtentContainsPoint(const kmAABBExtent* a, const kmVec3* p)
{
  if(p->x < a->point.x - a->extent.x || p->x > a->point.x + a->extent.x) return KM_FALSE;
  if(p->y < a->point.y - a->extent.y || p->y > a->point.y + a->extent.y) return KM_FALSE;
  if(p->z < a->point.z - a->extent.z || p->z > a->point.z + a->extent.z) return KM_FALSE;
  return KM_TRUE;
}

/* calculate closest squared length from segments */
kmScalar kmClosestPointFromSegments(const kmVec3 *a1, const kmVec3 *b1, const kmVec3 *a2, const kmVec3 *b2,
      kmScalar *s, kmScalar *t, kmVec3 *c1, kmVec3 *c2)
{
   kmVec3 d1, d2, r, cd;
   kmScalar a, e, f, c, b, denom;
   assert(s && t && c1 && c2);

   kmVec3Subtract(&d1, b1, a1);
   kmVec3Subtract(&d2, b2, a2);
   kmVec3Subtract(&r, a1, a2);
   a = kmVec3Dot(&d1, &d1);
   e = kmVec3Dot(&d2, &d2);
   f = kmVec3Dot(&d2, &r);

   /* check if either or both segments degenerate into points */
   if (a <= kmEpsilon && e <= kmEpsilon) {
      /* both segments degenerate into points */
      *s = *t = 0.0f;
      kmVec3Assign(c1, a1);
      kmVec3Assign(c2, a2);
      kmVec3Subtract(&cd, c1, c2);
      return kmVec3Dot(&cd, &cd);
   }
   if (a <= kmEpsilon) {
      /* first segment degenerates into a point */
      *s = 0.0f;
      *t = f / e; // s = 0 => t = (b*s + f) / e = f / e
      *t = kmClamp(*t, 0.0f, 1.0f);
   } else {
      c = kmVec3Dot(&d1, &r);
      if (e <= kmEpsilon) {
         /* second segment degenerates into a point */
         *t = 0.0f;
         *s = kmClamp(-c / a, 0.0f, 1.0f); // t = 0 => s = (b*t - c) / a = -c / a
      } else {
         /* the general nondegenerate case starts here */
         b = kmVec3Dot(&d1, &d2);
         denom = a*e-b*b; /* always nonnegative */

         /* if segments not parallel, compute closest point on L1 to L2, and
          * clamp to segment S1. Else pick arbitrary s (here 0) */
         if (denom != 0.0f) {
            *s = kmClamp((b*f - c*e) / denom, 0.0f, 1.0f);
         } else *s = 0.0f;

         /* compute point on L2 closest to S1(s) using
          * t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f) / e */
         *t = (b*(*s) + f) / e;

         /* if t in [0,1] done. Else clamp t, recompute s for the new value
          * of t using s = Dot((P2+D2*t)-P1,D1) / Dot(D1,D1)= (t*b - c) / a
          * and clamp s to [0, 1] */
         if (*t < 0.0f) {
            *t = 0.0f;
            *s = kmClamp(-c / a, 0.0f, 1.0f);
         } else if (*t > 1.0f) {
            *t = 1.0f;
            *s = kmClamp((b - c) / a, 0.0f, 1.0f);
         }
      }
   }

   kmVec3Add(c1, a1, &d1);
   kmVec3Scale(c1, c1, *s);
   kmVec3Add(c2, a2, &d2);
   kmVec3Scale(c2, c2, *t);
   kmVec3Subtract(&cd, c1, c2);
   return kmVec3Dot(&cd, &cd);
}

static kmScalar kmSqDistPointAABB(const kmVec3 *p, const kmAABB *aabb)
{
   kmScalar sqDist = 0.0f;
   if (p->x < aabb->min.x) sqDist += (aabb->min.x - p->x) * (aabb->min.x - p->x);
   if (p->x > aabb->max.x) sqDist += (p->x - aabb->max.x) * (p->x - aabb->max.x);
   if (p->y < aabb->min.y) sqDist += (aabb->min.y - p->y) * (aabb->min.y - p->y);
   if (p->y > aabb->max.y) sqDist += (p->y - aabb->max.y) * (p->y - aabb->max.y);
   if (p->z < aabb->min.z) sqDist += (aabb->min.z - p->z) * (aabb->min.z - p->z);
   if (p->z > aabb->max.z) sqDist += (p->z - aabb->max.z) * (p->z - aabb->max.z);
   return sqDist;
}

static kmScalar kmSqDistPointAABBExtent(const kmVec3 *p, const kmAABBExtent *aabbe)
{
   kmAABB aabb;
   aabb.min.x = aabbe->point.x - aabbe->extent.x;
   aabb.max.x = aabbe->point.x + aabbe->extent.x;
   aabb.min.y = aabbe->point.y - aabbe->extent.y;
   aabb.max.y = aabbe->point.y + aabbe->extent.y;
   aabb.min.z = aabbe->point.z - aabbe->extent.z;
   aabb.max.z = aabbe->point.z + aabbe->extent.z;
   return kmSqDistPointAABB(p, &aabb);
}

static kmBool kmAABBExtentIntersectsLine(const kmAABBExtent* a, const kmVec3* p1, const kmVec3* p2)
{
   /* d = (p2 - p1) * 0.5 */
   kmVec3 d;
   kmVec3Subtract(&d, p2, p1);
   kmVec3Scale(&d, &d, 0.5f);

   /* c = p1 + d - (min + max) * 0.5; */
   kmVec3 c;
   kmVec3 min;
   kmVec3Scale(&min, &a->extent, -1.0f);
   kmVec3Add(&c, &min, &a->extent);
   kmVec3Scale(&c, &c, 0.5f);
   kmVec3Subtract(&c, &d, &c);
   kmVec3Add(&c, &c, p1);

   /* ad = abs(d) */
   kmVec3 ad;
   kmVec3Abs(&ad, &d);

   /* alias for clarity */
   const kmVec3* e = &a->extent;

   if (fabs(c.x) > e->x + ad.x) return KM_FALSE;
   if (fabs(c.y) > e->y + ad.y) return KM_FALSE;
   if (fabs(c.z) > e->z + ad.z) return KM_FALSE;
   if (fabs(d.y * c.z - d.z * c.y) > e->y * ad.z + e->z * ad.y + kmEpsilon) return KM_FALSE;
   if (fabs(d.z * c.x - d.x * c.z) > e->z * ad.x + e->x * ad.z + kmEpsilon) return KM_FALSE;
   if (fabs(d.x * c.y - d.y * c.x) > e->x * ad.y + e->y * ad.x + kmEpsilon) return KM_FALSE;

   return KM_TRUE;
}

static kmBool kmAABBIntersectsLine(const kmAABB* a, const kmVec3* p1, const kmVec3* p2)
{
   kmAABBExtent aabbe;
   kmAABBToAABBExtent(&aabbe, a);
   return kmAABBExtentIntersectsLine(&aabbe, p1, p2);
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

static void kmAABBApplyVelocity(kmAABB *aabb, const kmVec3 *velocity)
{
   kmVec3Add(&aabb->min, &aabb->min, velocity);
   kmVec3Add(&aabb->max, &aabb->max, velocity);
}

static void kmAABBEApplyVelocity(kmAABBExtent *aabbe, const kmVec3 *velocity)
{
   kmVec3Add(&aabbe->point, &aabbe->point, velocity);
}

static void kmOBBApplyVelocity(kmOBB *obb, const kmVec3 *velocity)
{
   kmAABBEApplyVelocity(&obb->aabb, velocity);
}

static void kmSphereApplyVelocity(kmSphere *sphere, const kmVec3 *velocity)
{
   kmVec3Add(&sphere->point, &sphere->point, velocity);
}

kmSphere* kmSphereFromAABB(kmSphere *sphere, const kmAABB *aabb)
{
   kmAABBCentre(aabb, &sphere->point);
   sphere->radius = kmAABBDiameterX(aabb);
   sphere->radius = max(sphere->radius, kmAABBDiameterY(aabb));
   sphere->radius = max(sphere->radius, kmAABBDiameterZ(aabb));
   return sphere;
}

kmBool kmSphereIntersectsAABBExtent(const kmSphere *a, const kmAABBExtent *b)
{
   if(kmAABBExtentContainsPoint(b, &a->point)) return KM_TRUE;
   kmScalar distance = kmSqDistPointAABBExtent(&a->point, b);
   return (distance <= a->radius * a->radius);
}

kmBool kmSphereIntersectsAABB(const kmSphere *a, const kmAABB *b)
{
   kmScalar distance = kmSqDistPointAABB(&a->point, b);
   return (distance <= a->radius * a->radius);
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

kmBool kmSphereIntersectsCapsule(const kmSphere *a, const kmCapsule *b)
{
   kmScalar distance, radiusSum;
   distance = kmVec3LengthSqSegment(&b->pointA, &b->pointB, &a->point);
   radiusSum = a->radius + b->radius;
   return (distance <= radiusSum * radiusSum);
}

kmBool kmCapsuleIntersectsCapsule(const kmCapsule *a, const kmCapsule *b)
{
   kmScalar s, t;
   kmVec3 c1, c2;
   kmScalar distance, radiusSum;
   distance = kmClosestPointFromSegments(&a->pointA, &a->pointB, &b->pointA, &b->pointB, &s, &s, &c1, &c2);
   radiusSum = a->radius + b->radius;
   return (distance <= radiusSum * radiusSum);
}

kmBool kmSphereIntersectsPlane(const kmSphere *a, const kmPlane *b)
{
   kmVec3 n = {b->a, b->b, b->c};
   kmScalar distance;
   distance = kmVec3Dot(&a->point, &n) - b->d;
   return (abs(distance) <= a->radius);
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
   for (i = 0; i < 3; ++i) for (j = 0; j < 3; ++j)
      mat.mat[i+j*3] = kmVec3Dot(&a->orientation[i], &b->orientation[j]);

   /* bring translations into a's coordinate frame */
   kmVec3Subtract(&tmp, &a->aabb.point, &b->aabb.point);
   translation.x = kmVec3Dot(&tmp, &a->orientation[0]);
   translation.y = kmVec3Dot(&tmp, &a->orientation[1]);
   translation.z = kmVec3Dot(&tmp, &a->orientation[2]);

   /* compute common subexpressions. add in and epsilon term to
    * counteract arithmetic errors when two edges are parallel and
    * their cross product is (near) null. */
   for (i = 0; i < 3; ++i) for (j = 0; j < 3; ++j)
      absMat.mat[i+j*3] = abs(mat.mat[i+j*3]) + kmEpsilon;

   /* test axes L = A0, L = A1, L = A2 */
   for (i = 0; i < 3; ++i) {
      ra = (i==0?a->aabb.extent.x:i==1?a->aabb.extent.y:a->aabb.extent.z);
      rb = b->aabb.extent.x * absMat.mat[i+0*3] + b->aabb.extent.y * absMat.mat[i+1*3] + b->aabb.extent.z * absMat.mat[i+2*3];
      if (abs((i==0?translation.x:i==1?translation.y:translation.z)) > ra + rb) return KM_FALSE;
   }

   /* test axes L = B0, L = B1, L = B2 */
   for (i = 0; i < 3; ++i) {
      ra = a->aabb.extent.x * absMat.mat[0+i*3] + a->aabb.extent.y * absMat.mat[1+i*3] + a->aabb.extent.z * absMat.mat[2+i*3];
      rb = (i==0?b->aabb.extent.x:i==1?b->aabb.extent.y:b->aabb.extent.z);
      if (abs(translation.x * mat.mat[0+i*3] + translation.y * mat.mat[1+i*3] + translation.z * mat.mat[2+i*3]) > ra + rb) return KM_FALSE;
   }

   /* test axis L = A0 x B0 */
   ra = a->aabb.extent.y * absMat.mat[2+0*3] + a->aabb.extent.z * absMat.mat[1+0*3];
   rb = b->aabb.extent.y * absMat.mat[0+2*3] + b->aabb.extent.z * absMat.mat[0+1*3];
   if (abs(translation.z * mat.mat[1+0*3] - translation.y * mat.mat[2+0*3]) > ra + rb) return KM_FALSE;

   /* test axis L = A0 x B1 */
   ra = a->aabb.extent.y * absMat.mat[2+1*3] + a->aabb.extent.z * absMat.mat[1+1*3];
   rb = b->aabb.extent.x * absMat.mat[0+2*3] + b->aabb.extent.z * absMat.mat[0+0*3];
   if (abs(translation.z * mat.mat[1+1*3] - translation.y * mat.mat[2+1*3]) > ra + rb) return KM_FALSE;

   /* test axis L = A0 x B2 */
   ra = a->aabb.extent.y * absMat.mat[2+2*3] + a->aabb.extent.z * absMat.mat[1+2*3];
   rb = b->aabb.extent.x * absMat.mat[0+1*3] + b->aabb.extent.y * absMat.mat[0+0*3];
   if (abs(translation.z * mat.mat[1+2*3] - translation.y * mat.mat[2+2*3]) > ra + rb) return KM_FALSE;

   /* test axis L = A1 x B0 */
   ra = a->aabb.extent.x * absMat.mat[2+0*3] + a->aabb.extent.z * absMat.mat[0+0*3];
   rb = b->aabb.extent.y * absMat.mat[1+2*3] + b->aabb.extent.z * absMat.mat[1+1*3];
   if (abs(translation.x * mat.mat[2+0*3] - translation.z * mat.mat[0+0*3]) > ra + rb) return KM_FALSE;

   /* test axis L = A1 x B1 */
   ra = a->aabb.extent.x * absMat.mat[2+1*3] + a->aabb.extent.z * absMat.mat[0+1*3];
   rb = b->aabb.extent.x * absMat.mat[1+2*3] + b->aabb.extent.z * absMat.mat[1+0*3];
   if (abs(translation.x * mat.mat[2+1*3] - translation.z * mat.mat[0+1*3]) > ra + rb) return KM_FALSE;

   /* test axis L = A1 x B2 */
   ra = a->aabb.extent.x * absMat.mat[2+2*3] + a->aabb.extent.z * absMat.mat[0+2*3];
   rb = b->aabb.extent.x * absMat.mat[1+1*3] + b->aabb.extent.y * absMat.mat[1+0*3];
   if (abs(translation.x * mat.mat[2+2*3] - translation.z * mat.mat[0+2*3]) > ra + rb) return KM_FALSE;

   /* test axis L = A2 x B0 */
   ra = a->aabb.extent.x * absMat.mat[1+0*3] + a->aabb.extent.y * absMat.mat[0+0*3];
   rb = b->aabb.extent.y * absMat.mat[2+2*3] + b->aabb.extent.z * absMat.mat[2+1*3];
   if (abs(translation.y * mat.mat[0+0*3] - translation.x * mat.mat[1+0*3]) > ra + rb) return KM_FALSE;

   /* test axis L = A2 x B1 */
   ra = a->aabb.extent.x * absMat.mat[1+1*3] + a->aabb.extent.y * absMat.mat[0+1*3];
   rb = b->aabb.extent.x * absMat.mat[2+2*3] + b->aabb.extent.z * absMat.mat[2+0*3];
   if (abs(translation.y * mat.mat[0+1*3] - translation.x * mat.mat[1+1*3]) > ra + rb) return KM_FALSE;

   /* test axis L = A2 x B2 */
   ra = a->aabb.extent.x * absMat.mat[1+2*3] + a->aabb.extent.y * absMat.mat[0+2*3];
   rb = b->aabb.extent.x * absMat.mat[2+1*3] + b->aabb.extent.y * absMat.mat[2+0*3];
   if (abs(translation.y * mat.mat[0+2*3] - translation.x * mat.mat[1+2*3]) > ra + rb) return KM_FALSE;

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

kmBool kmAABBIntersectsSphere(const kmAABB *a, const kmSphere *b)
{
   return kmSphereIntersectsAABB(b, a);
}

kmBool kmAABBExtentIntersectsAABBExtent(const kmAABBExtent *a, const kmAABBExtent *b)
{
   if (abs(a->point.x - b->point.x) > (a->extent.x + b->extent.x)) return KM_FALSE;
   if (abs(a->point.y - b->point.y) > (a->extent.y + b->extent.y)) return KM_FALSE;
   if (abs(a->point.z - b->point.z) > (a->extent.z + b->extent.z)) return KM_FALSE;
   return KM_TRUE;
}

kmBool kmAABBExtentIntersectsSphere(const kmAABBExtent *a, const kmSphere *b)
{
   return kmSphereIntersectsAABBExtent(b, a);
}

kmBool kmAABBIntersectsAABBExtent(const kmAABB *a, const kmAABBExtent *b)
{
   if (a->min.x > b->point.x + b->extent.x || a->max.x < b->point.x - b->extent.x) return KM_FALSE;
   if (a->min.y > b->point.y + b->extent.y || a->max.y < b->point.y - b->extent.y) return KM_FALSE;
   if (a->min.z > b->point.z + b->extent.z || a->max.z < b->point.z - b->extent.z) return KM_FALSE;
   return KM_TRUE;
}

kmBool kmAABBExtentIntersectsAABB(const kmAABBExtent *a, const kmAABB *b)
{
   return kmAABBIntersectsAABBExtent(b, a);
}

kmBool kmAABBExtentIntersectsOBB(const kmAABBExtent *a, const kmOBB *b)
{
   kmOBB obb = {*a, {{1,0,0},{0,1,0},{0,0,1}}};
   return kmOBBIntersectsOBB(&obb, b);
}

kmBool kmOBBIntersectsAABBExtent(const kmOBB *a, const kmAABBExtent *b)
{
   return kmAABBExtentIntersectsOBB(b, a);
}

kmBool kmAABBIntersectsOBB(const kmAABB *a, const kmOBB *b)
{
   kmAABBExtent aabbExtent;
   kmAABBToAABBExtent(&aabbExtent, a);
   return kmAABBExtentIntersectsOBB(&aabbExtent, b);
}

kmBool kmOBBIntersectsAABB(const kmOBB *a, const kmAABB *b)
{
   return kmAABBIntersectsOBB(b, a);
}

kmBool kmAABBExtentIntersectsCapsule(const kmAABBExtent *a, const kmCapsule *b)
{
   /* Quick rejection test using the smallest (quickly calculatable) capsule the AABB fits inside*/
   kmCapsule smallestContainingCapsule = {{0, 0, 0}, {0, 0, 0}, 0};

   if(a->extent.x >= a->extent.y && a->extent.x >= a->extent.z)
   {
      smallestContainingCapsule.radius = sqrt(a->extent.y * a->extent.y + a->extent.z * a->extent.z);
      smallestContainingCapsule.pointA.x = a->point.x - a->extent.x;
      smallestContainingCapsule.pointB.x = a->point.x + a->extent.x;
   }
   else if(a->extent.y >= a->extent.x && a->extent.y >= a->extent.z)
   {
      smallestContainingCapsule.radius = sqrt(a->extent.x * a->extent.x + a->extent.z * a->extent.z);
      smallestContainingCapsule.pointA.y = a->point.y - a->extent.y;
      smallestContainingCapsule.pointB.y = a->point.y + a->extent.y;
   }
   else
   {
      smallestContainingCapsule.radius = sqrt(a->extent.x * a->extent.x + a->extent.y * a->extent.y);
      smallestContainingCapsule.pointA.z = a->point.z - a->extent.z;
      smallestContainingCapsule.pointB.z = a->point.z + a->extent.z;
   }

   if(!kmCapsuleIntersectsCapsule(&smallestContainingCapsule, b)) return KM_FALSE;

   /* Quick acceptance test for capsule line */
   if(kmAABBExtentIntersectsLine(a, &b->pointA, &b->pointB)) return KM_TRUE;

   /* Quick acceptance tests for capsule end spheres */
   kmSphere spa = {b->pointA, b->radius};
   if(kmAABBExtentIntersectsSphere(a, &spa)) return KM_TRUE;

   kmSphere spb = {b->pointB, b->radius};
   if(kmAABBExtentIntersectsSphere(a, &spb)) return KM_TRUE;



   return KM_TRUE;
}

/***
 * Collision manager
 ***/

/* internal */
typedef enum _CollisionType {
   COLLISION_AABB,
   COLLISION_AABBE,
   COLLISION_OBB,
   COLLISION_SPHERE,
   COLLISION_ELLIPSE,
   COLLISION_CAPSULE,
   COLLISION_LAST,
} _CollisionType;

typedef struct _CollisionPrimitive {
   /* shape union pointer */
   union {
      kmAABB *aabb;
      kmAABBExtent *aabbe;
      kmOBB *obb;
      kmSphere *sphere;
      kmEllipse *ellipse;
      kmCapsule *capsule;
      void *any;
   } shape;

   /* also a linked list */
   struct _CollisionPrimitive *next;

   /* user data */
   void *userData;

   /* type of shape */
   _CollisionType type;
} _CollisionPrimitive;

typedef struct _CollisionWorld {
   _CollisionPrimitive *primitives;
   void *userData;
} _CollisionWorld;

/* public */
typedef _CollisionWorld CollisionWorld;
typedef _CollisionPrimitive CollisionPrimitive;

struct CollisionInData;
struct CollisionOutData;
typedef void (*collisionResponseFunction)(const struct CollisionOutData *collision);
typedef int (*collisionTestFunction)(const struct CollisionInData *data, const CollisionPrimitive *collider);

/* public */
typedef struct CollisionInData {
   kmVec3 velocity;
   collisionResponseFunction response;
   collisionTestFunction test;
   const void *userData;
} CollisionInData;

typedef struct CollisionOutData {
   CollisionWorld *world;
   const CollisionPrimitive *collider;
   const kmVec3 *velocity;
   const kmVec3 *reflection;
   const void *userData;
} CollisionOutData;

typedef struct CollisionResponseInData {
   kmVec3 gravity;
} CollisionReponseInData;

typedef struct CollisionResponseOutData {
   kmVec3 position;
   char falling;
} CollisionResponseOutData;

/* internal */
typedef struct _CollisionPacket {
   /* shape union pointer */
   union {
      const kmAABB *aabb;
      const kmAABBExtent *aabbe;
      const kmOBB *obb;
      const kmSphere *sphere;
      const kmEllipse *ellipse;
      const kmCapsule *capsule;
      const void *any;
   } shape;

   /* input data we use to calculate out data */
   const CollisionInData *inData;

   /* internal data */
   kmScalar nearestDistance;

   /* type of shape */
   _CollisionType type;
} _CollisionPacket;

const char* _collisionTypeString(_CollisionType type)
{
   switch (type) {
      case COLLISION_AABB:return "AABB";
      case COLLISION_AABBE:return "AABBE";
      case COLLISION_OBB:return "OBB";
      case COLLISION_SPHERE:return "SPHERE";
      case COLLISION_ELLIPSE:return "ELLIPSE";
      case COLLISION_CAPSULE:return "CAPSULE";
      default:break;
   }
   return "INVALID TYPE";
}

static void _collisionPacketPrimitiveReflection(const _CollisionPacket *packet, const _CollisionPrimitive *primitive, void *shape, void *velocityFunction, void *testFunction, kmVec3 *outReflection)
{
   int i;
   kmVec3 reflection, inverseReflection;
   typedef kmBool (*_collisionTestFunc)(const void *a, const void *b);
   typedef void (*_velocityApplyFunc)(const void *a, const kmVec3 *velocity);
   _collisionTestFunc test = (_collisionTestFunc)testFunction;
   _velocityApplyFunc velocity = (_velocityApplyFunc)velocityFunction;

   /* state before collision */
   kmVec3Scale(&inverseReflection, &packet->inData->velocity, -1.0f);
   velocity(shape, &inverseReflection);

   /* find reflection for collision */
   memset(&reflection, 0, sizeof(kmVec3));
   for (i = 0; i < 2; ++i) {
      if (i == 0) reflection.y = packet->inData->velocity.y;
      else if (i == 1) reflection.x = packet->inData->velocity.x;
      else if (i == 2) reflection.z = packet->inData->velocity.z;

      kmVec3Scale(&inverseReflection, &reflection, -1.0f);
      velocity(shape, &reflection);

      if (test(shape, primitive->shape.any)) {
         if (i == 0) reflection.y *= -1.0f;
         else if (i == 1) reflection.x *= -1.0f;
         else if (i == 2) reflection.z *= -1.0f;
      } else {
         if (i == 0) reflection.y = 0.0f;
         else if (i == 1) reflection.x = 0.0f;
         else if (i == 2) reflection.z = 0.0f;
      }

      velocity(shape, &inverseReflection);
   }

   memcpy(outReflection, &reflection, sizeof(kmVec3));
}

static void _collisionPrimitiveFindPacketReflection(const _CollisionPrimitive *primitive, const _CollisionPacket *packet, void *testFunction, kmVec3 *outReflection)
{
   kmAABB aabb;
   kmAABBExtent aabbe;
   kmOBB obb;
   kmSphere sphere;
   memset(outReflection, 0, sizeof(kmVec3));

   switch (packet->type) {
      case COLLISION_AABB:
         memcpy(&aabb, packet->shape.any, sizeof(kmAABB));
         _collisionPacketPrimitiveReflection(packet, primitive, &aabb, kmAABBApplyVelocity, testFunction, outReflection);
         break;
      case COLLISION_AABBE:
         memcpy(&aabbe, packet->shape.any, sizeof(kmAABBExtent));
         _collisionPacketPrimitiveReflection(packet, primitive, &aabbe, kmAABBEApplyVelocity, testFunction, outReflection);
         break;
      case COLLISION_OBB:
         memcpy(&obb, packet->shape.any, sizeof(kmOBB));
         _collisionPacketPrimitiveReflection(packet, primitive, &obb, kmOBBApplyVelocity, testFunction, outReflection);
         break;
      case COLLISION_SPHERE:
         memcpy(&sphere, packet->shape.any, sizeof(kmSphere));
         _collisionPacketPrimitiveReflection(packet, primitive, &sphere, kmSphereApplyVelocity, testFunction, outReflection);
         break;
      case COLLISION_ELLIPSE:
      case COLLISION_CAPSULE:
      default:break;
   }
}

static void _collisionWorldPrimitiveTestWithPacket(CollisionWorld *object, const _CollisionPrimitive *primitive, _CollisionPacket *packet, void *testFunction)
{
   typedef kmBool (*_collisionTestFunc)(const void *a, const void *b);
   _collisionTestFunc test = (_collisionTestFunc)testFunction;

   if (!test(packet->shape.any, primitive->shape.any))
      return;

   // printf("%s COLLIDE %s\n", _collisionTypeString(packet->type), _collisionTypeString(primitive->type));

   kmVec3 reflection;
   _collisionPrimitiveFindPacketReflection(primitive, packet, testFunction, &reflection);

   CollisionOutData outData;
   outData.world = object;
   outData.collider = primitive;
   outData.velocity = &packet->inData->velocity;
   outData.reflection = &reflection;
   outData.userData = packet->inData->userData;
   if (packet->inData->response) packet->inData->response(&outData);
}

static void _collisionWorldCollideWithPacket(CollisionWorld *object, _CollisionPacket *packet)
{
   _CollisionPrimitive *p;

   void *testFunction[COLLISION_LAST][COLLISION_LAST];
   memset(testFunction, 0, sizeof(testFunction));

   /* aabb vs. x */
   testFunction[COLLISION_AABB][COLLISION_AABB] = kmAABBIntersectsAABB;
   testFunction[COLLISION_AABB][COLLISION_SPHERE] = kmAABBIntersectsSphere;
   testFunction[COLLISION_AABB][COLLISION_AABBE] = kmAABBIntersectsAABBExtent;
   testFunction[COLLISION_AABB][COLLISION_OBB] = kmAABBIntersectsOBB;

   /* aabbe vs. x */
   testFunction[COLLISION_AABBE][COLLISION_AABBE] = kmAABBExtentIntersectsAABBExtent;
   testFunction[COLLISION_AABBE][COLLISION_SPHERE] = kmAABBExtentIntersectsSphere;
   testFunction[COLLISION_AABBE][COLLISION_AABB] = kmAABBExtentIntersectsAABB;
   testFunction[COLLISION_AABBE][COLLISION_OBB] = kmAABBExtentIntersectsOBB;

   /* sphere vs. x */
   testFunction[COLLISION_SPHERE][COLLISION_SPHERE] = kmSphereIntersectsSphere;
   testFunction[COLLISION_SPHERE][COLLISION_AABB] = kmSphereIntersectsAABB;
   testFunction[COLLISION_SPHERE][COLLISION_AABBE] = kmSphereIntersectsAABBExtent;

   /* obb vs. x */
   testFunction[COLLISION_OBB][COLLISION_OBB] = kmOBBIntersectsOBB;
   testFunction[COLLISION_OBB][COLLISION_AABB] = kmOBBIntersectsAABB;
   testFunction[COLLISION_OBB][COLLISION_AABBE] = kmOBBIntersectsAABBExtent;

   for (p = object->primitives; p; p = p->next) {
      /* ask user if we should even bother testing */
      if (!packet->inData->test(packet->inData, p))
         continue;

      /* do the test */
      if (testFunction[packet->type][p->type]) {
         _collisionWorldPrimitiveTestWithPacket(object, p, packet, testFunction[packet->type][p->type]);
      } else {
         printf("-!- Test function between primitives not implemented!\n");
      }
   }
}

static void _collisionPrimitiveFree(_CollisionPrimitive *primitive)
{
   assert(primitive);
   switch (primitive->type) {
      default:
         IFDO(free, primitive->shape.any);
         break;
   }
   IFDO(free, primitive);
}

static _CollisionPrimitive* _collisionWorldAddPrimitive(CollisionWorld *object, _CollisionType type, void *shape, void *userData)
{
   _CollisionPrimitive *primitive, *p;
   assert(object && shape);

   if (!(primitive = calloc(1, sizeof(_CollisionPrimitive))))
      goto fail;

   primitive->type = type;
   primitive->shape.any = shape;
   primitive->userData = userData;

   if (!(p = object->primitives)) {
      object->primitives = primitive;
   } else {
      for (; p && p->next; p = p->next);
      p->next = primitive;
   }

   return primitive;

fail:
   IFDO(free, primitive);
   return NULL;
}

CollisionWorld* collisionWorldNew(void *userData)
{
   CollisionWorld *object;

   if (!(object = calloc(1, sizeof(CollisionWorld))))
      goto fail;

   object->userData = userData;
   return object;

fail:
   return NULL;
}

void collisionWorldFree(CollisionWorld *object)
{
   _CollisionPrimitive *p, *pn;
   assert(object);

   for (p = object->primitives; p; p = pn) {
      pn = p->next;
      _collisionPrimitiveFree(p);
   }

   free(object);
}

void* collisionWorldGetUserData(const CollisionWorld *object)
{
   assert(object);
   return object->userData;
}

CollisionPrimitive* collisionWorldAddEllipse(CollisionWorld *object, const kmEllipse *ellipse, void *userData)
{
   kmEllipse *ellipseCopy = NULL;
   _CollisionPrimitive *primitive = NULL;
   assert(object && ellipse);

   if (!(ellipseCopy = malloc(sizeof(kmEllipse))))
      goto fail;

   memcpy(ellipseCopy, ellipse, sizeof(kmEllipse));
   if (!(primitive = _collisionWorldAddPrimitive(object, COLLISION_ELLIPSE, ellipseCopy, userData)))
      goto fail;
   return primitive;

fail:
   IFDO(free, primitive);
   IFDO(free, ellipseCopy);
   return NULL;
}

CollisionPrimitive* collisionWorldAddAABB(CollisionWorld *object, const kmAABB *aabb, void *userData)
{
   kmAABB *aabbCopy = NULL;
   _CollisionPrimitive *primitive = NULL;
   assert(object && aabb);

   if (!(aabbCopy = malloc(sizeof(kmAABB))))
      goto fail;

   memcpy(aabbCopy, aabb, sizeof(kmAABB));
   if (!(primitive = _collisionWorldAddPrimitive(object, COLLISION_AABB, aabbCopy, userData)))
      goto fail;
   return primitive;

fail:
   IFDO(free, aabbCopy);
   return NULL;
}

CollisionPrimitive* collisionWorldAddAABBExtent(CollisionWorld *object, const kmAABBExtent *aabbe, void *userData)
{
   kmAABBExtent *aabbeCopy = NULL;
   _CollisionPrimitive *primitive = NULL;
   assert(object && aabbe);

   if (!(aabbeCopy = malloc(sizeof(kmAABBExtent))))
      goto fail;

   memcpy(aabbeCopy, aabbe, sizeof(kmAABBExtent));
   if (!(primitive = _collisionWorldAddPrimitive(object, COLLISION_AABBE, aabbeCopy, userData)))
      goto fail;
   return primitive;

fail:
   IFDO(free, aabbeCopy);
   return NULL;
}

CollisionPrimitive* collisionWorldAddSphere(CollisionWorld *object, const kmSphere *sphere, void *userData)
{
   kmSphere *sphereCopy = NULL;
   _CollisionPrimitive *primitive = NULL;
   assert(object && sphere);

   if (!(sphereCopy = malloc(sizeof(kmSphere))))
      goto fail;

   memcpy(sphereCopy, sphere, sizeof(kmSphere));
   if (!(primitive = _collisionWorldAddPrimitive(object, COLLISION_SPHERE, sphereCopy, userData)))
      goto fail;
   return primitive;

fail:
   IFDO(free, sphereCopy);
   return NULL;
}

void collisionWorldRemovePrimitive(CollisionWorld *object, CollisionPrimitive *primitive)
{
   _CollisionPrimitive *p;
   assert(object && primitive);

   if (primitive == (p = object->primitives)) {
      object->primitives = primitive->next;
   } else {
      for (; p && p->next != primitive; p = p->next);
      if (p) p->next = primitive->next;
      else object->primitives = NULL;
   }

   _collisionPrimitiveFree(primitive);
}

void collisionWorldCollideAABB(CollisionWorld *object, const kmAABB *aabb, const CollisionInData *data)
{
   static kmVec3 zero = {0,0,0};
   _CollisionPacket packet;
   if (kmVec3AreEqual(&data->velocity, &zero))
      return;

   packet.type = COLLISION_AABB;
   packet.inData = data;
   // packet.nearestDistance = FLT_MAX;
   packet.shape.aabb = aabb;
   _collisionWorldCollideWithPacket(object, &packet);
}

void collisionWorldCollideAABBExtent(CollisionWorld *object, const kmAABBExtent *aabbe, const CollisionInData *data)
{
   static kmVec3 zero = {0,0,0};
   _CollisionPacket packet;
   if (kmVec3AreEqual(&data->velocity, &zero))
      return;

   packet.type = COLLISION_AABBE;
   packet.inData = data;
   // packet.nearestDistance = FLT_MAX;
   packet.shape.aabbe = aabbe;
   _collisionWorldCollideWithPacket(object, &packet);
}

void collisionWorldCollideOBB(CollisionWorld *object, const kmOBB *obb, const CollisionInData *data)
{
   static kmVec3 zero = {0,0,0};
   _CollisionPacket packet;
   if (kmVec3AreEqual(&data->velocity, &zero))
      return;

   packet.type = COLLISION_OBB;
   packet.inData = data;
   // packet.nearestDistance = FLT_MAX;
   packet.shape.obb = obb;
   _collisionWorldCollideWithPacket(object, &packet);
}

void collisionWorldCollideSphere(CollisionWorld *object, const kmSphere *sphere, const CollisionInData *data)
{
   static kmVec3 zero = {0,0,0};
   _CollisionPacket packet;
   if (kmVec3AreEqual(&data->velocity, &zero))
      return;

   packet.type = COLLISION_SPHERE;
   packet.inData = data;
   // packet.nearestDistance = FLT_MAX;
   packet.shape.sphere = sphere;
   _collisionWorldCollideWithPacket(object, &packet);
}

void* collisionPrimitiveGetUserData(const CollisionPrimitive *primitive)
{
   assert(primitive);
   return primitive->userData;
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

// Decomposes kmMat3 to Euler kmVec3
// X == Roll, Y == Yaw, Z == Pitch
static kmVec3* kmMat3ToEuler(kmVec3 *pOut, kmMat3 *pIn)
{
   // Check for gimbal locks
   if (kmAlmostEqual(pIn->mat[2], -1.0f)) {
      pOut->x = 0.0f;
      pOut->y = kmPI*0.5f;
      pOut->z = pOut->x + atan2f(pIn->mat[3], pIn->mat[6]);
   } else if (kmAlmostEqual(pIn->mat[2], 1.0f)) {
      pOut->x = 0.0f;
      pOut->y = -kmPI*0.5f;
      pOut->z = -pOut->x + atan2f(-pIn->mat[3], -pIn->mat[6]);
   } else {
      float y1 = -asinf(pIn->mat[2]);
      float y2 = kmPI - y1;

      float x1 = atan2f(pIn->mat[5] / cos(y1), pIn->mat[8] / cos(y1));
      float x2 = atan2f(pIn->mat[5] / cos(y2), pIn->mat[8] / cos(y2));

      float z2 = atan2f(pIn->mat[1] / cos(y2), pIn->mat[0] / cos(y2));
      float z1 = atan2f(pIn->mat[1] / cos(y1), pIn->mat[0] / cos(y1));

      // Find out shortest rotation
      if ((abs(x1) + abs(y1) + abs(z1)) <= (abs(x2) + abs(y2) + abs(z2))) {
         pOut->x = x1;
         pOut->y = y1;
         pOut->z = z1;
      } else {
         pOut->x = x2;
         pOut->y = y2;
         pOut->z = z2;
      }
   }

   pOut->x = kmRadiansToDegrees(pOut->x);
   pOut->y = kmRadiansToDegrees(pOut->y);
   pOut->z = kmRadiansToDegrees(pOut->z);
   return pOut;
}

glhckObject* glhckCubeFromKazmathOBB(const kmOBB *obb)
{
   kmMat3 mat;
   kmVec3 rot;
   glhckObject *o = glhckCubeFromKazmathAABBExtent(&obb->aabb);
   kmOBBGetMat3(obb, &mat);
   kmMat3ToEuler(&rot, &mat);
   glhckObjectRotation(o, &rot);
   return o;
}

glhckObject* glhckSphereFromKazmathSphere(const kmSphere *sphere)
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

typedef enum PrimitiveType {
   PRIMITIVE_AABB,
   PRIMITIVE_AABBE,
   PRIMITIVE_OBB,
   PRIMITIVE_SPHERE,
   PRIMITIVE_LAST,
} PrimitiveType;

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

typedef struct ResponseTest {
   PrimitiveType type;
   union {
      kmAABB *aabb;
      kmAABBExtent *aabbe;
      kmOBB *obb;
      kmSphere *sphere;
      void *any;
   } shape;
   collisionResponseFunction response;
} ResponseTest;

typedef struct PrimObjPair {
   void *primitive;
   glhckObject *object;
} PrimObjPair;

static int shouldTestAgainst(const CollisionInData *data, const CollisionPrimitive *collider)
{
   void *ourShape = ((PrimObjPair*)data->userData)->primitive;
   void *colliderShape = collisionPrimitiveGetUserData(collider);
   return (colliderShape != ourShape);
}

static void aabbResponse(const CollisionOutData *collision)
{
   glhckObject *object = ((PrimObjPair*)collision->userData)->object;
   kmAABB *aabb = ((PrimObjPair*)collision->userData)->primitive;

   if (object) {
      glhckMaterial *mat = glhckObjectGetMaterial(object);
      if (mat) glhckMaterialDiffuseb(mat, 255, 0, 0, 255);
   }

   kmAABBApplyVelocity(aabb, collision->reflection);
}

static void aabbeResponse(const CollisionOutData *collision)
{
   ResponseTest *tests = (ResponseTest*)collisionWorldGetUserData(collision->world);
   glhckObject *object = ((PrimObjPair*)collision->userData)->object;
   kmAABBExtent *aabbe = ((PrimObjPair*)collision->userData)->primitive;

   if (object) {
      glhckMaterial *mat = glhckObjectGetMaterial(object);
      if (mat) glhckMaterialDiffuseb(mat, 255, 0, 0, 255);
   }

   if (collisionPrimitiveGetUserData(collision->collider) == tests[0].shape.any) {
      kmVec3 velocity;
      kmVec3Scale(&velocity, collision->reflection, -1);
      kmAABBApplyVelocity(tests[0].shape.aabb, &velocity);

      CollisionInData colData;
      memset(&colData, 0, sizeof(colData));
      memcpy(&colData.velocity, &velocity, sizeof(kmVec3));
      PrimObjPair pair;
      pair.primitive = tests[0].shape.any;
      pair.object = NULL;
      colData.response = tests[0].response;
      colData.test = shouldTestAgainst;
      colData.userData = &pair;
      collisionWorldCollideAABB(collision->world, tests[0].shape.aabb, &colData);
   }

   if (collisionPrimitiveGetUserData(collision->collider) == tests[2].shape.any) {
      kmVec3 velocity;
      kmVec3Scale(&velocity, collision->reflection, -1);
      kmSphereApplyVelocity(tests[2].shape.sphere, &velocity);

      CollisionInData colData;
      memset(&colData, 0, sizeof(colData));
      memcpy(&colData.velocity, &velocity, sizeof(kmVec3));
      PrimObjPair pair;
      pair.primitive = tests[2].shape.any;
      pair.object = NULL;
      colData.response = tests[2].response;
      colData.test = shouldTestAgainst;
      colData.userData = &pair;
      collisionWorldCollideSphere(collision->world, tests[2].shape.sphere, &colData);
   }

   kmAABBEApplyVelocity(aabbe, collision->reflection);
}

static void obbResponse(const CollisionOutData *collision)
{
   kmOBB *obb = ((PrimObjPair*)collision->userData)->primitive;
   ((PrimObjPair*)collision->userData)->primitive = &obb->aabb;
   aabbeResponse(collision);
}

static void sphereResponse(const CollisionOutData *collision)
{
   glhckObject *object = ((PrimObjPair*)collision->userData)->object;
   kmSphere *sphere = ((PrimObjPair*)collision->userData)->primitive;

   if (object) {
      glhckMaterial *mat = glhckObjectGetMaterial(object);
      if (mat) glhckMaterialDiffuseb(mat, 255, 0, 0, 255);
   }

   kmSphereApplyVelocity(sphere, collision->reflection);
}

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
   PrimitiveType primitive = PRIMITIVE_AABB;
   unsigned int i, textSize, font, testsInPrimitive = 0;
   char keyLock = 0, expectIntersection = 0, didIntersect = 0;
   char *primitiveStr = NULL, basicIntersections = 1;
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
            .orientation[0] = {0.71,-0.71,0},
            .orientation[1] = {0.71, 0.71,0},
            .orientation[2] = {0, 0,      1}}),
       .b = (&(kmOBB){.aabb = *aabbeTests[0].b,
            .orientation[0] = {0.71,-0.71,0},
            .orientation[1] = {0.71, 0.71,0},
            .orientation[2] = {0, 0,      1}}),
      .intersection = 0},
      {.a = (&(kmOBB){.aabb = *aabbeTests[1].a,
            .orientation[0] = {0.71,-0.71,0},
            .orientation[1] = {0.71, 0.71,0},
            .orientation[2] = {0, 0,      1}}),
       .b = (&(kmOBB){.aabb = *aabbeTests[1].b,
            .orientation[0] = {0.71,-0.71,0},
            .orientation[1] = {0.71, 0.71,0},
            .orientation[2] = {0, 0,      1}}),
      .intersection = 0},
      {.a = (&(kmOBB){.aabb = *aabbeTests[2].a,
            .orientation[0] = {0.71,-0.71,0},
            .orientation[1] = {0.71, 0.71,0},
            .orientation[2] = {0, 0,      1}}),
       .b = (&(kmOBB){.aabb = *aabbeTests[2].b,
            .orientation[0] = {0.71,-0.71,0},
            .orientation[1] = {0.71, 0.71,0},
            .orientation[2] = {0, 0,      1}}),
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

   ResponseTest responseTests[] = {
      {.type = PRIMITIVE_AABB,
       .shape.aabb = (&(kmAABB){.min = {0,0,0}, .max = {0,0,0}}),
       .response = aabbResponse},
      {.type = PRIMITIVE_AABBE,
       .shape.aabbe = (&(kmAABBExtent){.point = {0,0,0}, .extent = {0,0,0}}),
       .response = aabbeResponse},
#if 0
      {.type = PRIMITIVE_OBB,
       .shape.obb = (&(kmOBB){.aabb = {.point = {0,0,0}, .extent = {0,0,0}},
            .orientation[0] = {0,0,0},
            .orientation[1] = {0,0,0},
            .orientation[2] = {0,0,0}}),
       .response = obbResponse},
#endif
      {.type = PRIMITIVE_SPHERE,
       .shape.sphere = (&(kmSphere){.point = {0,0,0}, .radius = 0}),
       .response = sphereResponse},
   };
   ResponseTest responseTestsInitial[] = {
      {.type = PRIMITIVE_AABB,
       .shape.aabb = (&(kmAABB){.min = {WIDTH*0.2-50,HEIGHT/2-50,0}, .max = {WIDTH*0.2+50,HEIGHT/2+50,0}}),
       .response = aabbResponse},
      {.type = PRIMITIVE_AABBE,
       .shape.aabbe = (&(kmAABBExtent){.point = {WIDTH*0.4,HEIGHT/2,0}, .extent = {50,50,0}}),
       .response = aabbeResponse},
#if 0
      {.type = PRIMITIVE_OBB,
       .shape.obb = (&(kmOBB){.aabb = {.point = {WIDTH*0.6,HEIGHT/2,0}, .extent = {50,50,0}},
            .orientation[0] = {0.71,-0.71,0},
            .orientation[1] = {0.71, 0.71,0},
            .orientation[2] = {0, 0,      1}}),
       .response = obbResponse},
#endif
      {.type = PRIMITIVE_SPHERE,
       .shape.sphere = (&(kmSphere){.point = {WIDTH*0.8,HEIGHT/2,0}, .radius = 50}),
       .response = sphereResponse},
   };

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
   createPrimitives(i, sphereTests, glhckSphereFromKazmathSphere);
   kmAABB groundAABB = { .min = { 0, HEIGHT-80, 0 }, .max = { WIDTH, HEIGHT, 0 } };
   kmAABB wall1AABB = { .min = { -10, 40, 0 }, .max = { 5, HEIGHT-80, 0 } };

   i = 0;
   while(RUNNING && glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS) {
      glfwPollEvents();

      if (basicIntersections) {
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
      } else {
         CollisionWorld *world = collisionWorldNew(responseTests);
         collisionWorldAddAABB(world, &groundAABB, NULL);
         collisionWorldAddAABB(world, &wall1AABB, NULL);
         for (i = 0; i != LENGTH(responseTests); ++i) {
            ResponseTest *test = &responseTests[i];
            switch (test->type) {
               case PRIMITIVE_AABB:
                  collisionWorldAddAABB(world, test->shape.aabb, test->shape.aabb);
                  break;
               case PRIMITIVE_AABBE:
                  collisionWorldAddAABBExtent(world, test->shape.aabbe, test->shape.aabbe);
                  break;
               case PRIMITIVE_SPHERE:
                  collisionWorldAddSphere(world, test->shape.sphere, test->shape.sphere);
                  break;
            }
         }

         for (i = 0; i != LENGTH(responseTests); ++i) {
            ResponseTest *test = &responseTests[i];
            glhckObject *object = NULL;

            CollisionInData colData;
            colData.velocity = (kmVec3){ 0, 0.1, 0 };
            if (test->type == PRIMITIVE_AABBE) {
               colData.velocity.x += glfwGetKey(window, GLFW_KEY_RIGHT) * 0.1;
               colData.velocity.x -= glfwGetKey(window, GLFW_KEY_LEFT) * 0.1;
               colData.velocity.y -= glfwGetKey(window, GLFW_KEY_UP) * 0.2;
               colData.velocity.y += glfwGetKey(window, GLFW_KEY_DOWN) * 0.2;
            }

            switch (test->type) {
               case PRIMITIVE_AABB:
                  object = glhckCubeFromKazmathAABB(test->shape.aabb);
                  kmAABBApplyVelocity(test->shape.aabb, &colData.velocity);
                  break;
               case PRIMITIVE_AABBE:
                  object = glhckCubeFromKazmathAABBExtent(test->shape.aabbe);
                  kmAABBEApplyVelocity(test->shape.aabbe, &colData.velocity);
                  break;
               case PRIMITIVE_OBB:
                  object = glhckCubeFromKazmathOBB(test->shape.obb);
                  kmOBBApplyVelocity(test->shape.obb, &colData.velocity);
                  break;
               case PRIMITIVE_SPHERE:
                  object = glhckSphereFromKazmathSphere(test->shape.sphere);
                  kmSphereApplyVelocity(test->shape.sphere, &colData.velocity);
                  break;
               default:break;
            }
            if (!object) continue;

            glhckObjectDrawAABB(object, 1);
            glhckMaterial *mat = glhckMaterialNew(NULL);
            glhckMaterialDiffuseb(mat, 0, 255, 0, 255);
            glhckObjectMaterial(object, mat);
            glhckMaterialFree(mat);

            PrimObjPair pair;
            pair.primitive = test->shape.any;
            pair.object = object;
            colData.response = test->response;
            colData.test = shouldTestAgainst;
            colData.userData = &pair;

            switch (test->type) {
               case PRIMITIVE_AABB:
                  collisionWorldCollideAABB(world, test->shape.aabb, &colData);
                  break;
               case PRIMITIVE_AABBE:
                  collisionWorldCollideAABBExtent(world, test->shape.aabbe, &colData);
                  break;
               case PRIMITIVE_OBB:
                  collisionWorldCollideOBB(world, test->shape.obb, &colData);
                  break;
               case PRIMITIVE_SPHERE:
                  collisionWorldCollideSphere(world, test->shape.sphere, &colData);
                  break;
               default:break;
            }

            glhckObjectRender(object);
            glhckObjectFree(object);
         }

         glhckObject *object = glhckCubeFromKazmathAABB(&groundAABB);
         glhckMaterial *mat = glhckMaterialNew(NULL);
         glhckMaterialDiffuseb(mat, 255, 255, 255, 255);
         glhckObjectMaterial(object, mat);
         glhckMaterialFree(mat);
         glhckObjectRender(object);
         glhckObjectFree(object);

         object = glhckCubeFromKazmathAABB(&wall1AABB);
         mat = glhckMaterialNew(NULL);
         glhckMaterialDiffuseb(mat, 255, 255, 255, 255);
         glhckObjectMaterial(object, mat);
         glhckMaterialFree(mat);
         glhckObjectRender(object);
         glhckObjectFree(object);

         collisionWorldFree(world);
      }

      /* next intersection test */
      if (glfwGetKey(window, GLFW_KEY_SPACE)) {
         if (!keyLock && ++i >= testsInPrimitive) {
            if (++primitive >= PRIMITIVE_LAST) primitive = PRIMITIVE_AABB;
            i = 0;
         }
         keyLock = 1;
      } else if (glfwGetKey(window, GLFW_KEY_TAB)) {
         if (!keyLock) basicIntersections = !basicIntersections;

         for (i = 0; i != LENGTH(responseTests); ++i) {
            ResponseTest *test = &responseTests[i];
            switch (test->type) {
               case PRIMITIVE_AABB:
                  memcpy(test->shape.any, responseTestsInitial[i].shape.any, sizeof(kmAABB));
                  break;
               case PRIMITIVE_AABBE:
                  memcpy(test->shape.any, responseTestsInitial[i].shape.any, sizeof(kmAABBExtent));
                  break;
               case PRIMITIVE_OBB:
                  memcpy(test->shape.any, responseTestsInitial[i].shape.any, sizeof(kmOBB));
                  break;
               case PRIMITIVE_SPHERE:
                  memcpy(test->shape.any, responseTestsInitial[i].shape.any, sizeof(kmSphere));
                  break;
               default:break;
            }
         }

         i = 0;
         keyLock = 1;
      } else {
         keyLock = 0;
      }

      /* draw hud */
      glhckTextStash(text, font, textSize, 0, textSize, "colhck - coltest.c", NULL);
      glhckTextStash(text, font ,textSize, 0, textSize*2, "press [tab] to switch between primitive/response tests", NULL);

      if (basicIntersections) {
         glhckTextStash(text, font, textSize, 0, textSize*3, "press [space] for next test case", NULL);
         glhckTextStash(text, font, textSize, 0, textSize*5, primitiveStr, NULL);
         if (expectIntersection) {
            glhckTextStash(text, font, textSize, 0, textSize*6, "expect intersection: YES", NULL);
         } else {
            glhckTextStash(text, font, textSize, 0, textSize*6, "expect intersection: NO", NULL);
         }
         glhckTextRender(text);
         glhckTextClear(text);

         if (didIntersect == expectIntersection) {
            glhckTextColorb(text, 0, 255, 0, 255);
            glhckTextStash(text, font, textSize, 0, textSize*7, "CORRECT", NULL);
         } else {
            glhckTextColorb(text, 255, 0, 0, 255);
            glhckTextStash(text, font, textSize, 0, textSize*7, "WRONG", NULL);
         }
         glhckTextRender(text);
         glhckTextClear(text);
         glhckTextColorb(text, 255, 255, 255, 255);
      } else {
         glhckTextRender(text);
         glhckTextClear(text);
      }

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
