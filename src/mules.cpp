
#include "data.h"

double interpolateAlpha(Data2D& data, int faceId, char direction){
}

void calculateHighOrderFlux(Data2D& data, int faceId){
    Face2D& curFace = data.faces[faceId];
    if (faceId < data.nhorizontalFaces){
        if (curFace.v > 0){
            curFace.alphaFlux = curFace.v[CORRECTED_2] * data.dt * interpolateAlpha(data, curFace.id, 'forward');
        } else {
            curFace.alphaFlux = curFace.v[CORRECTED_2] * data.dt * interpolateAlpha(data, curFace.id, 'backward');
        }
    } else {
        if (curFace.u > 0){
            curFace.alphaFlux = curFace.u[CORRECTED_2] * data.dt * interpolateAlpha(data, curFace.id, 'forward');
        } else {
            curFace.alphaFlux = curFace.u[CORRECTED_2] * data.dt * interpolateAlpha(data, curFace.id, 'backward');
        }
    }
}


    