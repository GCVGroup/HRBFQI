#ifndef AXISKDTREE
#define AXISKDTREE

#include "PointSet.h"
#include <stdio.h>
#include <iostream>

#define SEARCH_MAX 10

class AxisKdTree;

class AxisKdCell{
public:
  
  AxisKdTree* _tree;
  int _start, _end;
  
  AxisKdCell* _fore;
  AxisKdCell* _back;
  
  char _s_axis;
  float _middle;
  
  bool _leaf;
  
  AxisKdCell(AxisKdTree* tree, int s, int e){
    _tree = tree;
    _start = s;
    _end = e;
    
    _leaf = true;
    
    _fore = NULL;
    _back = NULL;
  }
  
  ~AxisKdCell(){
    if(!_leaf){
      delete _fore;
      delete _back;
    }
  }
  
  bool split();
  void collectPointIndexInSphere(float c[3], float r);
  void countPointInBox(int& count, float min[3], float max[3]);
  void collectPointInBox(float (*point)[3], float *value,
                         int& count, float min[3], float max[3]);
  void getPointBound(float minP[3], float maxP[3], float min[3], float max[3]);
  
  void build(){
    if(split()){
      if(_fore != NULL){
        _fore->build();
      }
      if(_back != NULL){
        _back->build();
      }
    }
  }
};

class AxisKdTree{
public:
  PointSet* _ps;
  AxisKdCell* _root;
  int *_index_list;
  int _listN;
  
  AxisKdTree(PointSet* ps){
    _ps = ps;
    
    int N = ps->point_N;
    _root = new AxisKdCell(this, 0, N);
    
    _index_list = new int[N];
    
    _root->build();
  }

  ~AxisKdTree(){
    delete _root;
    delete[] _index_list;
  }
  
  void collectPointIndexInSphere(int* &index_list, int &listN, float c[3], float r){
    _listN = 0;
    _root->collectPointIndexInSphere(c, r);
    
    index_list = _index_list;
    listN = _listN;
  }
  
  int countPointInBox(float min[3], float max[3]){
    int count = 0;
    _root->countPointInBox(count, min, max);
    return count;
  }
  
  void collectPointInBox(float (*point)[3], float *value, float min[3], float max[3]){
    int count = 0;
    _root->collectPointInBox(point, value, count, min, max);
  }
  
  void getPointBound(float min[3], float max[3]){
    float minP[3] = {max[0], max[1], max[2]};
    float maxP[3] = {min[0], min[1], min[2]};
    _root->getPointBound(minP, maxP, min, max);
    
    min[0] = minP[0];
    min[1] = minP[1];
    min[2] = minP[2];
    
    max[0] = maxP[0];
    max[1] = maxP[1];
    max[2] = maxP[2];
  }
};

#endif
