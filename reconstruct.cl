#define LEFT -1
#define RIGHT 1
#define DOWN -1
#define UP 1

void ReconstructX1(const __local REAL* restrict primTile, 
                   const int iTile, const int jTile,
                   REAL* restrict primEdge,
                   const int dir)
{
    REAL left, mid, right, slope;

    for (int var=0; var<DOF; var++) {
        left = primTile[INDEX_LOCAL(iTile-1, jTile, var)];
        mid = primTile[INDEX_LOCAL(iTile, jTile, var)];
        right = primTile[INDEX_LOCAL(iTile+1, jTile, var)];
        
        slope = SlopeLim(left, mid, right);

        primEdge[var] = mid + dir*0.5*slope;
    }
}

void ReconstructX2(const __local REAL* restrict primTile, 
                   const int iTile, const int jTile,
                   REAL* restrict primEdge,
                   const int dir)
{
    REAL left, mid, right, slope;

    for (int var=0; var<DOF; var++) {
        left = primTile[INDEX_LOCAL(iTile, jTile-1, var)];
        mid = primTile[INDEX_LOCAL(iTile, jTile, var)];
        right = primTile[INDEX_LOCAL(iTile, jTile+1, var)];
        
        slope = SlopeLim(left, mid, right);

        primEdge[var] = mid + dir*0.5*slope;
    }

}
