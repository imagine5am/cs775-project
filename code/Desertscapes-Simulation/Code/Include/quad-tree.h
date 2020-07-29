#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <vec.h>
#include <vector>


class Quad;

/*
    Holds cell data. Quad constructor to populate most fields in the starting.
    Later these can be updated.
*/

//struct floatvals{
//    float sediments;
//    float bedrock;
//    float vegetation;
//    floatvals() {
//        sediments = 0.0;
//        bedrock = 0.0;
//        vegetation =0.0;
//    }
//};

class CellData
{

public:
    Vector2 botLeft;        //bounds
    Vector2 topRight;       //bounds
    Vector2 centre;

    float sizex;
    float sizey;

    std::map<char, float> vals;

    Quad* qd;


    CellData()
    {
        botLeft = Vector2(-1, -1);
        topRight = Vector2(-1, -1);
        centre = Vector2(-1, -1);
        sizex = -1.0;
        sizey = -1.0;

        qd = NULL;
       /* vals.insert({'b',0.0});
        vals.insert({ 's',0.0 });
        vals.insert({ 'v',0.0 });*/
       
        vals['b'] = 0.0; //bedrock
        vals['s'] = 0.0; //sediment
        vals['v'] = 0.0; //vegetation
        //std::cerr << vals.size() << std::endl;
    }
};

/*
- The following is the quad tree structure.
- It should contain mostly the algorithmic implementation.
- With the exception of a one time write process in the constructor,
  simulation specifications must go in cd.
- Bounds of the quad to be fetched from cd.
-Assume full 4-ary tree : All children or no children
*/
class Quad
{
    //Vector2 type bounds are needed. - Shivam

    // Contains details of node 
    

    // Children of this tree 



public:
    static std::vector<std::vector<Quad*>>* qroots;
    static int nx, ny;
    CellData* cd;
    Quad* tLChild;
    Quad* tRChild;
    Quad* bLChild;
    Quad* bRChild;
    Quad* parent;

    int level;
    Vector2 bLBound;
    Vector2 tRBound;

    //Quad()
    //{
    //    parent = NULL;
    //    tLChild = NULL;
    //    tRChild = NULL;
    //    bLChild = NULL;
    //    bRChild = NULL;
    //    level = 0;
    //}

 

    Quad(Vector2 botL, Vector2 topR, Quad* par, bool isLeaf,int lvl)
    {
        //std::cerr << "chkpt4.cons 1" << std::endl;
        //qroots = new std::vector<std::vector<Quad*>>();
        parent = par;
        tLChild = NULL;
        tRChild = NULL;
        bLChild = NULL;
        bRChild = NULL;

        bLBound = botL;
        tRBound = topR;

        level = lvl;

        //isLeaf = true;
        if (isLeaf) {
            cd = new CellData();


            cd->topRight = topR;
            cd->botLeft = botL;
            cd->sizex = (topR - botL).x;





            cd->sizey = (topR - botL).y;
            cd->centre = (topR + botL) / 2;
        }
        //else cd = NULL;

        
        //std::cerr << "chkpt4.cons end" << std::endl;
    }

    CellData* fetchData();

    CellData* search(Vector2 pt);                 
    std::vector<CellData*> leftNeighbors(float top, float bot);
    std::vector<CellData*> leftNeighbors();
    std::vector<CellData*> rightNeighbors(float top, float bot);
    std::vector<CellData*> rightNeighbors();
    std::vector<CellData*> topNeighbors(float top, float bot);
    std::vector<CellData*> topNeighbors();
    std::vector<CellData*> bottomNeighbors(float top, float bot);
    std::vector<CellData*> bottomNeighbors();

    std::vector<CellData*> allNeighbors();


    bool withinBounds(Vector2);                 //Shivam - not mandatory

    bool istLChild();
    bool istRChild();
    bool isbLChild();
    bool isbRChild();
    bool isTChild();
    bool isBChild();
    bool isLChild();
    bool isRChild();

    std::vector<CellData*> allBottomCellData();
    std::vector<CellData*> allBottomCellData(float top, float bot);
    std::vector<CellData*> allTopCellData();
    std::vector<CellData*> allTopCellData(float top, float bot);
    std::vector<CellData*> allLeftCellData();
    std::vector<CellData*> allLeftCellData(float top, float bot);
    std::vector<CellData*> allRightCellData();
    std::vector<CellData*> allRightCellData(float top, float bot);
    

    std::vector<CellData*> allLeafCellData();
    bool isLeaf();

    void insert( int numLevels, std::map<char, float> vals);
};



inline void Quad::insert( int numLevels, std::map<char, float> vals) {
    //max lvl is number of levels - 1

    Vector2 botL = bLBound;
    Vector2 topR = tRBound;
    Vector2 mid = (botL + topR) / 2;


    if (level == numLevels - 2) {

        //Quad temp(Vector2(botL.x, mid.y), Vector2(mid.x, topR.y), this, true, level + 1);
        
        // = Quad(Vector2(botL.x, mid.y), Vector2(mid.x, topR.y), this, true, level + 1);
        //this->tLChild = &Quad( Vector2(botL.x,mid.y), Vector2(mid.x, topR.y), this, true, level + 1);
        
        //this->tLChild = &temp;

        tLChild = new Quad(Vector2(botL.x, mid.y), Vector2(mid.x, topR.y), this, true, level + 1);

        tLChild->cd->vals = vals;

        tLChild->cd->qd = tLChild;

        tRChild = new Quad(mid, topR, this, true, level + 1);
        tRChild->cd->vals = vals;
        tRChild->cd->qd = tRChild;

        bLChild = new Quad(botL, mid, this, true, level + 1);
        bLChild->cd->vals = vals;
        bLChild->cd->qd = bLChild;

        bRChild = new Quad(Vector2(mid.x, botL.y), Vector2(topR.x, mid.y), this, true, level + 1);
        bRChild->cd->vals = vals;
        bRChild->cd->qd = bRChild;
    }

    if (level < numLevels - 2) {

        tLChild = new Quad(Vector2(botL.x, mid.y), Vector2(mid.x, topR.y), this, false, level + 1);
        tLChild->insert(numLevels,vals);

        tRChild = new Quad(mid, topR, this, false, level + 1);
        tRChild->insert(numLevels,vals);

        bLChild = new Quad(botL, mid, this, false, level + 1);
        bLChild->insert(numLevels,vals);

        bRChild = new Quad(Vector2(mid.x, botL.y), Vector2(topR.x, mid.y), this, false, level + 1);
        bRChild->insert(numLevels,vals);
    
    }


}


//Assume leaf quad is searching
inline std::vector<CellData*> Quad::leftNeighbors(float top, float bot) {

   //Assertion: We only enter this function after confirming the neighbors are in same tree.


    if(this->parent->tRChild == this)
        return this->parent->tLChild->allRightCellData(top, bot);

    if (this->parent->bRChild == this)
        return this->parent->bLChild->allRightCellData(top, bot);

    //bot-up checks
    return this->parent->leftNeighbors(top, bot);


}

inline std::vector<CellData*> Quad::leftNeighbors() {

    float top = tRBound.y;
    float bot = bLBound.y;
    std::vector<CellData*> retVal;
    if (float(int(bLBound.x)) == bLBound.x) {  //adjacent root

        if (int(bLBound.x) > 0){
            return qroots->at(int(bLBound.x) - 1).at(int(bLBound.y) )->allRightCellData(top, bot);

        }
        else {
            return retVal; //no neighnor :(
        }
    }


    return leftNeighbors(top,bot);
}

//inline std::vector<CellData*> Quad::leftNeighbors() {
//    std::vector<CellData*> retVal;
//    if (this->level == 0) {
//        //This is the root. Find the left roots.
//        // Shivam Todo - get left_root of this.
//        Quad* left_root;
//        std::vector<CellData*> leafs = left_root->allRightCellData();
//        retVal.insert(leafs.end(), leafs.begin(), leafs.end());
//    }
//    else {
//        if (istRChild()) {
//            std::vector<CellData*> leafs = this->parent->tLChild->allRightCellData();
//            retVal.insert(leafs.end(), leafs.begin(), leafs.end());
//        }
//        else if (isbRChild()) {
//            std::vector<CellData*> leafs = this->parent->bLChild->allRightCellData();
//            retVal.insert(leafs.end(), leafs.begin(), leafs.end());
//        }
//        else if (istLChild()) {
//            Quad* left_node;
//            if (this->parent->level == 0) {
//                left_node = NULL;                       // Get left root
//                if (left_node->isLeaf()) {
//                    retVal.push_back(left_node->cd);
//                }
//                else {
//                    std::vector<CellData*> leafs = left_node->tRChild->allRightCellData();
//                    retVal.insert(leafs.end(), leafs.begin(), leafs.end());
//                }
//            }
//            else {
//                Quad* super_parent = this->parent->parent;
//                if (super_parent->istRChild()) {
//                    std::vector<CellData*> leafs = super_parent->tLChild->allRightCellData();
//                    retVal.insert(leafs.end(), leafs.begin(), leafs.end());
//                }
//                else if (super_parent->isbRChild()) {
//                    std::vector<CellData*> leafs = super_parent->bLChild->allRightCellData();
//                    retVal.insert(leafs.end(), leafs.begin(), leafs.end());
//                }
//                else if (super_parent->istLChild()) {
//
//                }
//                else if (super_parent->istLChild()) {
//
//                }
//            }
//        }
//        else if (isbLChild()) {
//
//        }
//        /*
//        if (istLChild()) {
//            if (this->parent->level == 0) {
//                // Shivam Todo - get left_root of this->parent
//                Quad* left_root;
//                if(left_root->tRChild != NULL) {
//                    std::vector<CellData*> leafs = left_root->tRChild->allRightCellData();
//                    retVal.insert(leafs.end(), leafs.begin(), leafs.end());
//                } else {
//                    retVal.push_back(left_root->cd);
//                }
//            } else {
//            }
//            // bottom-right most value of top-left neighbor of this->parent
//            // right most values in TR cell of left neighbor of this->parent
//            // top-right most value in BR cell of left neighbor of this->parent
//        } else if (istRChild()) {
//            std::vector<CellData*> leafs = this->parent->tLChild->allRightCellData();
//            retVal.insert(leafs.end(), leafs.begin(), leafs.end());
//            // bottom-right most value of top-left neighbor of this->parent
//            // right values of this->parent->tLChild;
//            // top-right most value of this->parent->bLChild;
//        } else if (isbLChild()) {
//        }
//        else if (isbRChild()) {
//        }
//        */
//    }
//    return retVal;
//}

inline std::vector<CellData*> Quad::rightNeighbors(float top, float bot) {
    //Assertion: We only enter this function after confirming the neighbors are in same tree.

    if (istLChild())
        return this->parent->tRChild->allLeftCellData(top, bot);
    else if (isbLChild())
        return this->parent->bRChild->allLeftCellData(top, bot);

    //bot-up checks
    return this->parent->rightNeighbors(top, bot);
}

inline std::vector<CellData*> Quad::rightNeighbors() {
 

    float top = tRBound.y;
    float bot = bLBound.y;
    std::vector<CellData*> retVal;
    if (float(int(tRBound.x)) == tRBound.x) {  //adjacent root
        if (int(tRBound.x) < nx)
            return qroots->at(int(bLBound.x) + 1).at(int(bLBound.y))->allLeftCellData(top, bot);
        else
            return retVal; //no neighnor :(
    }
    return rightNeighbors(top, bot);
}

inline std::vector<CellData*> Quad::topNeighbors(float top, float bot) {
    //Assertion: We only enter this function after confirming the neighbors are in same tree.

    if (isbLChild())
        return this->parent->tLChild->allBottomCellData(top, bot);
    else if (isbRChild())
        return this->parent->tRChild->allBottomCellData(top, bot);

    //bot-up checks
    return this->parent->topNeighbors(top, bot);
}

inline std::vector<CellData*> Quad::topNeighbors() {


    float top = tRBound.x;
    float bot = bLBound.x;
    std::vector<CellData*> retVal;
    if (float(int(tRBound.y)) == tRBound.y) {  //adjacent root
        if (int(tRBound.y) < ny)
            return qroots->at(int(bLBound.x) ).at(int(bLBound.y) + 1)->allBottomCellData(top, bot);
        else
            return retVal; //no neighnor :(
    }
    return topNeighbors(top, bot);
}

inline std::vector<CellData*> Quad::bottomNeighbors(float top, float bot) {
    //Assertion: We only enter this function after confirming the neighbors are in same tree.

    if (istLChild())
        return this->parent->bLChild->allTopCellData(top, bot);
    else if (istRChild())
        return this->parent->bRChild->allTopCellData(top, bot);

    //bot-up checks
    return this->parent->bottomNeighbors(top, bot);
}

inline std::vector<CellData*> Quad::bottomNeighbors() {

    float top = tRBound.x;
    float bot = bLBound.x;
    std::vector<CellData*> retVal;
    if (float(int(bLBound.y)) == bLBound.y) {  //adjacent root
        if (int(bLBound.y) > 0)
            return qroots->at(int(bLBound.x)).at(int(bLBound.y) - 1)->allTopCellData(top, bot);
        else
            return retVal; //no neighnor :(
    }
    return bottomNeighbors(top, bot);
}
inline bool Quad::istLChild() {
    if (this->parent->tLChild == this)
        return true;
    return false;
}
inline bool Quad::istRChild() {
    if (this->parent->tRChild == this)
        return true;
    return false;
}
inline bool Quad::isbLChild() {
    if (this->parent->bLChild == this)
        return true;
    return false;
}
inline bool Quad::isbRChild() {
    if (this->parent->bRChild == this)
        return true;
    return false;
}
inline bool Quad::isTChild() {
    return istLChild() || istRChild();
}
inline bool Quad::isBChild() {
    return isbLChild() || isbRChild();
}
inline bool Quad::isLChild() {
    return istLChild() || isbLChild();
}
inline bool Quad::isRChild() {
    return istRChild() || isbRChild();
}

inline std::vector<CellData*> Quad::allBottomCellData(float top, float bot) {
    if (tRBound.x <= top && bLBound.x >= bot) {
        return allBottomCellData();
    }

    std::vector<CellData*> retVal;
    if (isLeaf()) {
        retVal.push_back(this->cd);
    }
    else {
        std::queue<Quad*> queue_obj;
        if (this->bRChild->bLBound.x < top) //<= : corners
            queue_obj.push(this->bRChild);

        if (this->bLChild->tRBound.x > bot)
            queue_obj.push(this->bLChild);

        while (!queue_obj.empty()) {
            Quad* elem = queue_obj.front();
            queue_obj.pop();
            if (elem->isLeaf()) {
                retVal.push_back(elem->cd);
            }
            else {
                if (elem->bRChild->bLBound.x < top)
                    queue_obj.push(elem->bRChild);

                if (elem->bLChild->tRBound.x > bot)
                    queue_obj.push(elem->bLChild);
            }
        }
    }
    return retVal;
}

inline std::vector<CellData*> Quad::allBottomCellData() {
    std::vector<CellData*> retVal;
    if (isLeaf()) {
        retVal.push_back(this->cd);
    }
    else {
        std::queue<Quad*> queue_obj;
        queue_obj.push(this->bLChild);
        queue_obj.push(this->bRChild);
        while (!queue_obj.empty()) {
            Quad* elem = queue_obj.front();
            queue_obj.pop();
            if (elem->isLeaf()) {
                retVal.push_back(elem->cd);
            }
            else {
                queue_obj.push(elem->bLChild);
                queue_obj.push(elem->bRChild);
            }
        }
    }
    return retVal;
}

inline std::vector<CellData*> Quad::allTopCellData(float top, float bot) {
    if (tRBound.x <= top && bLBound.x >= bot) {
        return allBottomCellData();
    }

    std::vector<CellData*> retVal;
    if (isLeaf()) {
        retVal.push_back(this->cd);
    }
    else {
        std::queue<Quad*> queue_obj;
        if (this->tRChild->bLBound.x < top) //<= : corners
            queue_obj.push(this->tRChild);

        if (this->tLChild->tRBound.x > bot)
            queue_obj.push(this->tLChild);

        while (!queue_obj.empty()) {
            Quad* elem = queue_obj.front();
            queue_obj.pop();
            if (elem->isLeaf()) {
                retVal.push_back(elem->cd);
            }
            else {
                if (elem->tRChild->bLBound.x < top)
                    queue_obj.push(elem->tRChild);

                if (elem->tLChild->tRBound.x > bot)
                    queue_obj.push(elem->tLChild);
            }
        }
    }
    return retVal;
}

inline std::vector<CellData*> Quad::allTopCellData() {
    std::vector<CellData*> retVal;
    if (isLeaf()) {
        retVal.push_back(this->cd);
    }
    else {
        std::queue<Quad*> queue_obj;
        queue_obj.push(this->tLChild);
        queue_obj.push(this->tRChild);
        while (!queue_obj.empty()) {
            Quad* elem = queue_obj.front();
            queue_obj.pop();
            if (elem->isLeaf()) {
                retVal.push_back(elem->cd);
            }
            else {
                queue_obj.push(elem->tLChild);
                queue_obj.push(elem->tRChild);
            }
        }
    }
    return retVal;
}

inline std::vector<CellData*> Quad::allRightCellData(float top, float bot) {

    if (tRBound.y <= top && bLBound.y >= bot) {
        return allRightCellData();
    }

    std::vector<CellData*> retVal;
    if (isLeaf()) {
        retVal.push_back(this->cd);
    }
    else {
        std::queue<Quad*> queue_obj;
        if (this->tRChild->bLBound.y < top) //<= : corners
            queue_obj.push(this->tRChild);

        if (this->bRChild->tRBound.y > bot)
            queue_obj.push(this->bRChild);

        while (!queue_obj.empty()) {
            Quad* elem = queue_obj.front();
            queue_obj.pop();
            if (elem->isLeaf()) {
                retVal.push_back(elem->cd);
            }
            else {
                if (elem->tRChild->bLBound.y < top)
                    queue_obj.push(elem->tRChild);

                if (elem->bRChild->tRBound.y > bot)
                    queue_obj.push(elem->bRChild);
            }
        }
    }
    return retVal;
}

inline std::vector<CellData*> Quad::allRightCellData() {
    std::vector<CellData*> retVal;
    if (isLeaf()) {
        retVal.push_back(this->cd);
    }
    else {
        std::queue<Quad*> queue_obj;
        queue_obj.push(this->tRChild);
        queue_obj.push(this->bRChild);
        while (!queue_obj.empty()) {
            Quad* elem = queue_obj.front();
            queue_obj.pop();
            if (elem->isLeaf()) {
                retVal.push_back(elem->cd);
            }
            else {
                queue_obj.push(elem->tRChild);
                queue_obj.push(elem->bRChild);
            }
        }
    }
    return retVal;
}

inline std::vector<CellData*> Quad::allLeftCellData(float top, float bot) {
    if (tRBound.y <= top && bLBound.y >= bot) {
        return allLeftCellData();
    }

    std::vector<CellData*> retVal;
    if (isLeaf()) {
        retVal.push_back(this->cd);
    }
    else {
        std::queue<Quad*> queue_obj;
        if (this->tLChild->bLBound.y < top) //<= : corners
            queue_obj.push(this->tLChild);

        if (this->bLChild->tRBound.y > bot)
            queue_obj.push(this->bLChild);

        while (!queue_obj.empty()) {
            Quad* elem = queue_obj.front();
            queue_obj.pop();
            if (elem->isLeaf()) {
                retVal.push_back(elem->cd);
            }
            else {
                if (elem->tLChild->bLBound.y < top)
                    queue_obj.push(elem->tLChild);

                if (elem->bLChild->tRBound.y > bot)
                    queue_obj.push(elem->bLChild);
            }
        }
    }
    return retVal;
}

inline std::vector<CellData*> Quad::allLeftCellData() {
    std::vector<CellData*> retVal;
    if (isLeaf()) {
        retVal.push_back(this->cd);
    }
    else {
        std::queue<Quad*> queue_obj;
        queue_obj.push(this->tLChild);
        queue_obj.push(this->bLChild);
        while (!queue_obj.empty()) {
            Quad* elem = queue_obj.front();
            queue_obj.pop();
            if (elem->isLeaf()) {
                retVal.push_back(elem->cd);
            }
            else {
                queue_obj.push(elem->tLChild);
                queue_obj.push(elem->bLChild);
            }
        }
    }
    return retVal;
}



inline std::vector<CellData*> Quad::allLeafCellData() {
    std::vector<CellData*> retVal;
    if (isLeaf()) {
        retVal.push_back(this->cd);
    }
    else {
        std::queue<Quad*> queue_obj;
        queue_obj.push(this->tLChild);
        queue_obj.push(this->tRChild);
        queue_obj.push(this->bLChild);
        queue_obj.push(this->bRChild);

        while (!queue_obj.empty()) {
            Quad* elem = queue_obj.front();
            queue_obj.pop();
            if (elem->isLeaf()) {
                retVal.push_back(elem->cd);
            }
            else {
                queue_obj.push(elem->tLChild);
                queue_obj.push(elem->tRChild);
                queue_obj.push(elem->bLChild);
                queue_obj.push(elem->bRChild);
            }
        }
    }
    return retVal;
}

inline bool Quad::isLeaf() {
    return (this->cd != NULL);// ? true : false;
}


inline std::vector<CellData*> Quad::allNeighbors() {
    std::vector<CellData*> ret, temp;
    ret = leftNeighbors();
    temp = rightNeighbors();
    ret.insert(ret.end(), temp.begin(), temp.end());
    
    temp = topNeighbors();
    ret.insert(ret.end(), temp.begin(), temp.end());
    
    temp = bottomNeighbors();
    ret.insert(ret.end(), temp.begin(), temp.end());

    return ret;
}

inline CellData* Quad::fetchData() {
    return cd;
}

/*pick any quad and search for a Vector2 point. Returns cell data*/
inline CellData* Quad::search(Vector2 pt) {

    //root check
    if ( (int(pt.x) != int(bLBound.x)) || (int(pt.y) != int(bLBound.y))  ) {
        return qroots->at(int(pt.x)).at(int(pt.y))->search(pt);
    }

    if(isLeaf())
        return cd;

    if (pt.x <= tRChild->bLBound.x) {        //Search Left
        if (pt.y <= tRChild->bLBound.y) {
            return bLChild->search(pt);
        }
        else return tLChild->search(pt);
    }
    else {                                  //Search Right
        if (pt.y <= tRChild->bLBound.y) {
            return bRChild->search(pt);
        }
        else return tRChild->search(pt);
    }

    return NULL;
}

