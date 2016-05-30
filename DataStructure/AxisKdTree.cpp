#include "AxisKdTree.h"

bool AxisKdCell::split(){
  if(_end-_start < SEARCH_MAX)
    return false;
  
  float min[3], max[3];
  _tree->_ps->getBound(min, max, _start, _end);
  if(max[0] - min[0] > max[1] - min[1]){
    _s_axis = 0;
    _middle = 0.5f*(max[0] + min[0]);
  }
  else{
    _s_axis = 1;
    _middle = 0.5f*(max[1] + min[1]);
  }
  if(max[2] - min[2] > max[_s_axis] - min[_s_axis]){
    _s_axis = 2;
    _middle = 0.5f*(max[2] + min[2]);
  }
  
  float (*point)[3] = _tree->_ps->point;
  
  int i = _start;
  int j = _end-1;
  while(i <= j){
    float p = point[i][_s_axis];
    if(p > _middle)
      i++;
    else{
      _tree->_ps->swapIndex(i, j);
      j--;
    }
  }
  
  if(_start == i || i == _end)
    return false;
  
  _leaf = false;
  _fore = new AxisKdCell(_tree, _start, i);
  _back = new AxisKdCell(_tree, i, _end);
  
  return true;
}

void AxisKdCell::collectPointIndexInSphere(float c[3], float r){
  if(_leaf){
     float (*point)[3] = _tree->_ps->point;
     double r2 = r*r;
     for(int i=_start; i<_end; i++){
       float *p = point[i];
       float vx = c[0] - p[0];
       float vy = c[1] - p[1];
       float vz = c[2] - p[2];
       if(vx*vx+vy*vy+vz*vz < r2){
         _tree->_index_list[_tree->_listN++] = i;
       }
     }
   }
   else{
     float d = c[_s_axis] - _middle;
     if(d+r >= 0)
       _fore->collectPointIndexInSphere(c, r);
     if(d-r <= 0)
       _back->collectPointIndexInSphere(c, r);
   }
}

void AxisKdCell::countPointInBox(int& count, float min[3], float max[3]){
  if(_leaf){
     float (*point)[3] = _tree->_ps->point;
     for(int i=_start; i<_end; i++){
       float *p = point[i];
       if(p[0] < min[0] || p[0] > max[0] ||
          p[1] < min[1] || p[1] > max[1] ||
          p[2] < min[2] || p[2] > max[2])
         continue;
       count++;
     }
   }
   else{
     if(_fore != NULL && max[_s_axis] >= _middle)
       _fore->countPointInBox(count, min, max);
     
     if(_back != NULL && min[_s_axis] <= _middle)
       _back->countPointInBox(count, min, max);
   }
}

//Is not used.
void AxisKdCell::collectPointInBox(float (*point)[3], float *value,
                                   int& count, float min[3], float max[3]){
  if(_leaf){
     for(int i=_start; i<_end; i++){
       float *p = _tree->_ps->point[i];
       if(p[0] < min[0] || p[0] > max[0] ||
          p[1] < min[1] || p[1] > max[1] ||
          p[2] < min[2] || p[2] > max[2])
         continue;
       point[count][0] = p[0];
       point[count][1] = p[1];
       point[count][2] = p[2];
       //value[count] = _tree->_ps->_value[i];
       count++;
     }
  }
   else{
     if(_fore != NULL && max[_s_axis] >= _middle)
       _fore->collectPointInBox(point, value, count, min, max);
     
     if(_back != NULL && min[_s_axis] <= _middle)
       _back->collectPointInBox(point, value, count, min, max);
   }
}

void AxisKdCell::getPointBound(float minP[3], float maxP[3], float min[3], float max[3]){
  if(_leaf){
     for(int i=_start; i<_end; i++){
       float *p = _tree->_ps->point[i];
       if(p[0] < min[0] || p[0] > max[0] ||
          p[1] < min[1] || p[1] > max[1] ||
          p[2] < min[2] || p[2] > max[2])
         continue;
       if(p[0] < minP[0])
         minP[0] = p[0];
       if(p[0] > maxP[0])
         maxP[0] = p[0];
       
       if(p[1] < minP[1])
         minP[1] = p[1];
       if(p[1] > maxP[1])
         maxP[1] = p[1];
       
       if(p[2] < minP[2])
         minP[2] = p[2];
       if(p[2] > maxP[2])
         maxP[2] = p[2];
     }
  }
   else{
     if(_fore != NULL && max[_s_axis] >= _middle)
       _fore->getPointBound(minP, maxP, min, max);
     
     if(_back != NULL && min[_s_axis] <= _middle)
       _back->getPointBound(minP, maxP, min, max);
   }
}
