/* *
 *  Created on: 15-11-2012
 *      Author: miguel
 *      Editado: Fernando Santolaya 03/10/2025
 */
#include "K2tree.h"
#include <libcdsBasics.h>
#include <limits>
#include <queue>
#include "../Util/Factory.h"

using namespace std;


//Operacioens con representación de árbol
//creación de la representación
/*
void initVectorCuadrantSizes(int tamMatrix) {
    auto heightA = log(tamMatrix) / log(K);
    auto n = (size_t) pow(K, heightA);
    sizes.resize(heightA + 1);
    for (auto i = 0; i <= heightA; i++) {
        sizes[i] = n / pow(K, i);
    }
}*/
L_NODE *createL_Node() {
    L_NODE *resp = (L_NODE *) malloc(sizeof(L_NODE));
    resp->data = 0;
    resp->child = NULL;
    return resp;
}

lkt *llenaK2tree(string ruta, int elevate) {
    lkt *tmp = createLKTree(elevate); //16 ERA ANTES
    string linea1;
    ifstream set1;
    set1.open(ruta);
    int z, i = 0;
    string token;
    int aux[2];
    while (getline(set1, linea1)) {
        //Ahora leemos la linea y obtenemos los puntos
        istringstream ss(linea1);
        z = 0;
        while (getline(ss, token, ',')) {
            aux[z] = atoi(token.c_str());
            z++;
        }
        insertNode(tmp, aux[0], aux[1], i);
        i++;
    }
    set1.close();

    return tmp;
}

lkt *createLKTree(uint maxlevels) {
    lkt *tree = (lkt *) malloc(sizeof(lkt));
    tree->root = createL_Node();
    tree->max_Level = maxlevels;
    tree->numberNodes = 0;
    tree->numberLeaves = 0;
    tree->numberTotalLeaves = 0;
    return tree;
}

//delete a sub-tree recursively
int destroyLNode(L_NODE *node) {
    int resp = 0;
    if (node == NULL) {
        return resp;
    } else {
        if (node->child != NULL) {
            //delete child recursively
            for (int i = 0; i < K * K; i++) {
                resp += destroyLNode(node->child[i]);
                node->child[i] = NULL;
            }
            free(node->child);
            node->child = NULL;
            resp += sizeof(L_NODE *) * K * K + sizeof(L_NODE *) * K * K;
        }
        free(node);
        node = NULL;
        resp += sizeof(L_NODE) + sizeof(L_NODE *);
    }
    return resp;
}

//elimina la reprecentación temporal de un lktree.
int destroyLKTree(lkt *tree) {
    int resp = 0;
    if (tree == NULL)
        return resp;
    //Elimino la raiz
    resp = destroyLNode(tree->root);
    tree->root = NULL;
    //elimino el nodo cabecera.
    free(tree);
    tree = NULL;
    return resp;
}

void insertNode(lkt *tree, int x, int y, uint label) {
    uint *pt = (uint *) malloc(sizeof(uint));
    *pt = label;
    _insertNode(tree, x, y, pt);
}

void _insertNode(lkt *tree, int x, int y, uint *labelArray) {
    uint i, node = 0;
    unsigned long long int div_level;
    int l = 0;
    L_NODE *n = tree->root;
    while (l <= tree->max_Level) {
        div_level = pow(K, tree->max_Level - l);
        node = (x / div_level) * K + y / div_level;
        if (l == tree->max_Level) {
            if (n->data == 0) {
                tree->numberLeaves++;
            }
            n->data = n->data | (0x1 << node); //1;
        } else {
            if (n->child == NULL) {
                if (l < tree->max_Level - 1) {
                    tree->numberNodes += K * K;
                } else {
                    tree->numberTotalLeaves += K * K;
                }
                n->child = (L_NODE **) calloc(K * K, sizeof(L_NODE *));
                for (i = 0; i < K * K; i++) {
                    n->child[i] = NULL;
                }
            }
            if (n->child[node] == NULL) {
                n->child[node] = createL_Node();
            }

            //n apunta al nodo insertado (intermedio)
            n = n->child[node];
        }
        x = x % div_level;
        y = y % div_level;
        l++;
    }
}

int saveMREP2ToFile(MREP2 *ktree, std::ofstream &ofs) {
    if (!ktree) return -1;

    // Guardar los campos del MREP2
    ofs.write(reinterpret_cast<const char *>(&ktree->maxLevel), sizeof(int));
    ofs.write(reinterpret_cast<const char *>(&ktree->numberOfNodes), sizeof(size_t));
    ofs.write(reinterpret_cast<const char *>(&ktree->numberOfEdges), sizeof(size_t));

    // Guardar el bitmap BT
    if (ktree->bt) {
        size_t btSize = ktree->bt->getSize();
        ofs.write(reinterpret_cast<const char *>(&btSize), sizeof(size_t)); // Guardar tamaño del bitmap
        ktree->bt->save(ofs); // Usar el método de guardado de la clase
    } else {
        size_t btSize = 0;
        ofs.write(reinterpret_cast<const char *>(&btSize), sizeof(size_t)); // Guardar tamaño cero si no existe
    }

    // Guardar el bitmap BN
    if (ktree->bn) {
        size_t bnSize = ktree->bn->getSize();
        ofs.write(reinterpret_cast<const char *>(&bnSize), sizeof(size_t)); // Guardar tamaño del bitmap
        ktree->bn->save(ofs); // Usar el método de guardado de la clase
    } else {
        size_t bnSize = 0;
        ofs.write(reinterpret_cast<const char *>(&bnSize), sizeof(size_t)); // Guardar tamaño cero si no existe
    }

    // Guardar el bitmap BL
    if (ktree->bl) {
        size_t blSize = ktree->bl->getSize();
        ofs.write(reinterpret_cast<const char *>(&blSize), sizeof(size_t)); // Guardar tamaño del bitmap
        ktree->bl->save(ofs); // Usar el método de guardado de la clase
    } else {
        size_t blSize = 0;
        ofs.write(reinterpret_cast<const char *>(&blSize), sizeof(size_t)); // Guardar tamaño cero si no existe
    }

    // Guardar la tabla de divisiones de nivel
    ofs.write(reinterpret_cast<const char *>(ktree->div_level_table),
              sizeof(unsigned long long int) * (ktree->maxLevel));

    return 0;
}

int loadMREP2FromFile(MREP2 **ktree, std::ifstream &ifs) {
    *ktree = (MREP2 *) malloc(sizeof(MREP2));
    if (!*ktree) {
        perror("Error al asignar memoria para MREP2");
        return -1;
    }

    // Leer los campos básicos
    if (!ifs.read(reinterpret_cast<char *>(&(*ktree)->maxLevel), sizeof(int))) {
        perror("Error al leer maxLevel");
        free(*ktree);
        return -1;
    }
    if (!ifs.read(reinterpret_cast<char *>(&(*ktree)->numberOfNodes), sizeof(size_t))) {
        perror("Error al leer numberOfNodes");
        free(*ktree);
        return -1;
    }
    if (!ifs.read(reinterpret_cast<char *>(&(*ktree)->numberOfEdges), sizeof(size_t))) {
        perror("Error al leer numberOfEdges");
        free(*ktree);
        return -1;
    }

    //creammos el arreglo con los tamaños de los lados
    //initVectorCuadrantSizes((*ktree)->numberOfNodes);

    // Cargar el bitmap BT
    size_t btSize;
    if (!ifs.read(reinterpret_cast<char *>(&btSize), sizeof(size_t))) {
        perror("Error al leer el tamaño del bitmap BT");
        free(*ktree);
        return -1;
    }
    if (btSize > 0) {
        (*ktree)->bt = BitSequence::load(ifs); // Asegúrate de que esta función maneje la memoria correctamente
    } else {
        (*ktree)->bt = nullptr;
    }

    // Cargar el bitmap BN
    size_t bnSize;
    if (!ifs.read(reinterpret_cast<char *>(&bnSize), sizeof(size_t))) {
        perror("Error al leer el tamaño del bitmap BN");
        free(*ktree);
        return -1;
    }
    if (bnSize > 0) {
        (*ktree)->bn = BitSequence::load(ifs);
    } else {
        (*ktree)->bn = nullptr;
    }

    // Cargar el bitmap BL
    size_t blSize;
    if (!ifs.read(reinterpret_cast<char *>(&blSize), sizeof(size_t))) {
        perror("Error al leer el tamaño del bitmap BL");
        free(*ktree);
        return -1;
    }
    if (blSize > 0) {
        (*ktree)->bl = BitSequence::load(ifs);
    } else {
        (*ktree)->bl = nullptr;
    }

    // Cargar la tabla de divisiones de nivel
    (*ktree)->div_level_table = (unsigned long long int *)
            malloc(sizeof(unsigned long long int) * ((*ktree)->maxLevel));
    if (!(*ktree)->div_level_table) {
        perror("Error al asignar memoria para div_level_table");
        free(*ktree);
        return -1;
    }
    if (!ifs.read(reinterpret_cast<char *>((*ktree)->div_level_table),
                  sizeof(unsigned long long int) * ((*ktree)->maxLevel))) {
        perror("Error al leer div_level_table");
        free((*ktree)->div_level_table);
        free(*ktree);
        return -1;
    }

    return 0;
}

MREP2 *createK2tree(string dataset1, int elevate, uint cantPuntos) {
    //NO destruye lkt, hay que hacerlo desde afuera.
    //Estrategia general:
    //usando el arbol se construye el ktree de los puntos, y se construye un
    //arbol para las etiquetas.  Luego se construye el ktree para las estiquetas
    //con dicho arbol.

    //en este caso los valores de Numero de objetos y numero de arcos no se
    //si son necesarios, hay que revisar.
    int tamMatrix = pow(2, elevate);
    lkt *tree = llenaK2tree(dataset1, elevate);
    MREP2 *ktree = (MREP2 *) malloc(sizeof(MREP2));
    ktree->maxLevel = tree->max_Level;
    ktree->numberOfNodes = tamMatrix;
    ktree->numberOfEdges = cantPuntos;
    size_t bits_BT_len = tree->numberNodes;
    size_t bits_BN_len = tree->numberTotalLeaves;
    size_t bits_LI_len = tree->numberLeaves * K * K;

    //initVectorCuadrantSizes(tamMatrix);

    uint *bits_BT = new uint[uint_len(bits_BT_len, 1)]();
    uint *bits_BN = new uint[uint_len(bits_BN_len, 1)]();
    uint *bits_LI = new uint[uint_len(bits_LI_len, 1)]();

    uint k, j, queuecont, conttmp, node, pos = 0;

    unsigned long long div_level;
    int i;
    char isroot = 1;
    //Queue *q = createEmptyQueue();
    std::queue<L_NODE *> q;
    L_NODE *subTree;
    q.push(tree->root);
    queuecont = 1;

    for (i = 0; i < tree->max_Level; i++) {
        //i:iteración por nivel hasta el penúltimo
        conttmp = 0;
        div_level = pow(K, tree->max_Level - i);
        for (k = 0; k < queuecont; k++) {
            //por elemento de la cola de un mismo nivel
            subTree = (L_NODE *) q.front();
            q.pop();
            if ((subTree != NULL) && (subTree->child != NULL)) {
                for (j = 0; j < K * K; j++) {
                    node = j;
                    conttmp++;
                    q.push(subTree->child[node]);
                }
                if (!isroot)
                    bitset(bits_BT, pos);
            }
            if (!isroot)
                pos++;
            isroot = 0;
        }
        queuecont = conttmp;
    }

    ktree->bt = PlainBitSequenceFactory(bits_BT, bits_BT_len);
    pos = 0;
    size_t pos_inf = 0;

    //	fprintf(stderr,"Empezando bitmap de hojas utiles\n");
    while (!q.empty()) {
        subTree = (L_NODE *) q.front();
        q.pop();
        if ((subTree != NULL) && (subTree->data)) {
            //si no es 0
            bitset(bits_BN, pos);
            //hay que separar las listas

            for (i = 0; i < K * K; i++) {
                if ((subTree->data) & (0x1 << i)) {
                    bitset(bits_LI, pos_inf);
                }
                pos_inf++;
            }
        }
        pos++;
    }
    // destroyQueue(q);
    ktree->bn = PlainBitSequenceFactory(bits_BN, bits_BN_len);
    ktree->bl = PlainBitSequenceFactory(bits_LI, bits_LI_len);

    ktree->div_level_table = (long long uint *) malloc(
        sizeof(long long uint) * (ktree->maxLevel + 1));
    for (i = 0; i <= ktree->maxLevel; i++)
        ktree->div_level_table[i] = pow(K, ktree->maxLevel - i);

    delete[] bits_BT;
    bits_BT = NULL;
    delete[] bits_BN;
    bits_BN = NULL;
    delete[] bits_LI;
    bits_LI = NULL;
    return ktree;
}

void destroyK2tree(MREP2 *rep) {
    //destroyRepresentation(rep->ktree);
    delete rep->bl;
    delete rep->bn;
    delete rep->bt;
    free(rep->div_level_table);
    rep->bl = NULL;
    rep->bn = NULL;
    rep->bt = NULL;
    rep->div_level_table = NULL;
    free(rep);
    rep = NULL;
}

MREP2 *loadK2treeFromFile(const char *filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        perror("Error abriendo el archivo para cargar");
        return nullptr;
    }

    MREP2 *mrep = (MREP2 *) malloc(sizeof(MREP2));
    if (!mrep) {
        perror("Error al asignar memoria para MREP2");
        return nullptr;
    }

    // Leer el MREP2
    if (loadMREP2FromFile(&mrep, ifs) != 0) {
        perror("Error al cargar MREP2 desde archivo");
        free(mrep);
        return nullptr;
    }

    return mrep;
}

int saveK2treeToFile(MREP2 *snapshot, const char *filename) {
    if (!snapshot || !filename) return -1;

    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        perror("Error abriendo el archivo para guardar el Snapshot");
        return -1;
    }
    // Guardar el MREP2 (árbol k2tree)
    if (saveMREP2ToFile(snapshot, ofs) != 0) {
        ofs.close();
        return -1;
    }
    ofs.close();
    return 0;
}

void firstLeaf(MREP2 *rep, uint p1, uint p2, uint q1, uint q2,
               uint dp, uint dq, int x, int l, uint &n, Point &p) {
    uint i = 0, j, leaf;
    uint y, p1new, p2new, q1new, q2new;
    unsigned long int divlevel;
    if (l == rep->maxLevel) {
        //recorrido por el bitarray leavesInf
        leaf = x + i * p1;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                leaf = x + j;
                if (rep->bl->access(leaf)) {
                    //Point p;
                    p.setX(dp + i);
                    p.setY(dq + j);
                    return;
                    //dp + i => posición X del objeto econtrado
                    //dq + j => posición Y del objeto encontrado.
                    //snap->labels->getObjects(rep->bl->rank1(leaf), dp + i,dq + j, Oid, X, Y, n);
                }
            }
            leaf += K;
        }
    }

    if ((l == rep->maxLevel - 1)
        && (rep->bn->access(x - rep->bt->getLength()))) {
        //recorrido por el bitarray bn
        y = (rep->bn->rank1(x - rep->bt->getLength()) - 1) * K * K;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                firstLeaf(rep, 0, 0, 0, 0, dp + i, dq + j,
                          y + K * i + j, l + 1, n, p);
            }
        }
    }
    if ((x == -1) || ((l < rep->maxLevel - 1) && (rep->bt->access(x)))) {
        //recorrido por el bitarray bt
        y = (x == -1) ? 0 : rep->bt->rank1(x) * K * K;
        divlevel = rep->div_level_table[l + 1];
        for (i = p1 / divlevel; i <= p2 / divlevel; i++) {
            p1new = 0;
            if (i == p1 / divlevel)
                p1new = p1 % divlevel;
            p2new = divlevel - 1;
            if (i == p2 / divlevel)
                p2new = p2 % divlevel;
            for (j = q1 / divlevel; j <= q2 / divlevel; j++) {
                q1new = 0;
                if (j == q1 / divlevel)
                    q1new = q1 % divlevel;
                q2new = divlevel - 1;
                if (j == q2 / divlevel)
                    q2new = q2 % divlevel;
                firstLeaf(rep, p1new, p2new, q1new, q2new,
                          dp + divlevel * i, dq + divlevel * j, y + K * i + j,
                          l + 1, n, p);
            }
        }
    }
}

 int isLeaf(MREP2 *mrep, long int index) {
    //index está mapeado como si la implementación fuera un único bitmap
    //con la finalidad de encapsular las diferentes implementaciones existentes
    //del k^2-tree
    //verdadero si es una hoja de último nivel sino falso.
    //si index es -1 entonces es la raíz
    return (index != -1)
           && ((mrep->bt->getLength() + mrep->bn->getLength()) <= index);
}

 int posInLeaf(MREP2 *mrep, long int index) {
    return index - (mrep->bt->getLength() + mrep->bn->getLength());
}

 bool esCero(MREP2 *mrep, long int x) {
    long int btLen = mrep->bt->getLength();
    long int bnLen = mrep->bn->getLength();
    if (x == -1) {
        return false;
    } else if (x < btLen) {
        return mrep->bt->access(x) == 0;
    } else if (x < (btLen + bnLen)) {
        return mrep->bn->access(x - btLen) == 0;
    } else {
        //es una hoja de último nivel.
        return mrep->bl->access(x - (btLen + bnLen)) == 0;
    }
}

 long int firstChild(MREP2 *mrep, long int x) {
    long int btLen = mrep->bt->getLength();
    long int bnLen = mrep->bn->getLength();
    long int btOnes = mrep->bt->countOnes();
    if (x == -1) {
        return 0;
    } else if (x < btLen) {
        return (mrep->bt->access(x) == 0) ? -1 : mrep->bt->rank1(x) * K * K;
    } else if (x < btLen + bnLen) {
        /*		cout << endl << "btOnes:" << btOnes << " x in bt:" << x - btLen
         << " rank1 de x en bt:" << mrep->bn->rank1(x - btLen) << endl;
         cout << ((btOnes + (mrep->bn->rank1(x - btLen))) * K * K) << endl;
         */
        return (mrep->bn->access(x - btLen) == 0)
                   ?
                   //error: -1 : (btLen + (mrep->bn->rank1(x - btLen) * K * K));
                   -1
                   : ((btOnes + (mrep->bn->rank1(x - btLen))) * K * K);
    } else {
        //es una hoja de último nivel.
        return -1;
    }
}

void recursiveGetRank2(MREP2 *snap, uint p1, uint p2, uint q1, uint q2,
                       uint dp, uint dq, int x, int l, uint *&Oid, uint *&X,
                       uint *&Y, uint &n) {
    MREP2 *rep = snap;
    uint i = 0, j, leaf;
    uint y, p1new, p2new, q1new, q2new;
    unsigned long int divlevel;
    if (l == rep->maxLevel) {
        //recorrido por el bitarray leavesInf
        leaf = x + i * p1;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                leaf = x + j;
                if (rep->bl->access(leaf)) {
                    //dp + i => posición X del objeto econtrado
                    //dq + j => posición Y del objeto encontrado.
                    // snap->labels->getObjects(rep->bl->rank1(leaf), dp + i,
                    //                        dq + j, Oid, X, Y, n);
                    X[n] = dp + i;
                    Y[n] = dq + j;
                    n++;
                }
            }
            leaf += K;
        }
    }

    if ((l == rep->maxLevel - 1)
        && (rep->bn->access(x - rep->bt->getLength()))) {
        //recorrido por el bitarray bn
        y = (rep->bn->rank1(x - rep->bt->getLength()) - 1) * K * K;
        for (i = p1; i <= p2; i++) {
            for (j = q1; j <= q2; j++) {
                recursiveGetRank2(snap, 0, 0, 0, 0, dp + i, dq + j,
                                  y + K * i + j, l + 1, Oid, X, Y, n);
            }
        }
    }
    if ((x == -1) || ((l < rep->maxLevel - 1) && (rep->bt->access(x)))) {
        //recorrido por el bitarray bt
        y = (x == -1) ? 0 : rep->bt->rank1(x) * K * K;
        divlevel = rep->div_level_table[l + 1];
        for (i = p1 / divlevel; i <= p2 / divlevel; i++) {
            p1new = 0;
            if (i == p1 / divlevel)
                p1new = p1 % divlevel;
            p2new = divlevel - 1;
            if (i == p2 / divlevel)
                p2new = p2 % divlevel;
            for (j = q1 / divlevel; j <= q2 / divlevel; j++) {
                q1new = 0;
                if (j == q1 / divlevel)
                    q1new = q1 % divlevel;
                q2new = divlevel - 1;
                if (j == q2 / divlevel)
                    q2new = q2 % divlevel;
                recursiveGetRank2(snap, p1new, p2new, q1new, q2new,
                                  dp + divlevel * i, dq + divlevel * j, y + K * i + j,
                                  l + 1, Oid, X, Y, n);
            }
        }
    }
}

void rangeQueryHDOriginal2(MREP2 *snap, uint p1, uint p2, uint q1, uint q2,
                           uint *&Oid, uint *&X, uint *&Y, uint &n) {
    Oid = new uint[snap->numberOfEdges]();
    X = new uint[snap->numberOfEdges]();
    Y = new uint[snap->numberOfEdges]();
    n = 0;
    recursiveGetRank2(snap, p1, p2, q1, q2, 0, 0, -1, -1, Oid, X, Y, n);
}

std::vector<Point> extractPointK2tree(MREP2 *J, int nelementos) {
    uint *Oid1;
    uint *X1;
    uint *Y1;
    uint m1 = pow(K, J->maxLevel + 1) - 1;
    uint n1;
    // Llama a la función rangeQueryHDOriginal para obtener los puntos
    rangeQueryHDOriginal2(J, 0, m1, 0, m1, Oid1, X1, Y1, n1);
    std::vector<Point> Ar(n1);
    // Llena el vector con los puntos obtenidos
    for (uint i = 0; i < n1; i++) {
        Point p1;
        p1.setX(X1[i]);
        p1.setY(Y1[i]);
        Ar[i] = p1;
    }
    // Retorna el vector de puntos
    return Ar;
}

size_t sizeMREP2(MREP2 *rep) {
    size_t totalByte = 0;
    if (rep != NULL) {
        totalByte = sizeof(MREP2);
        totalByte += sizeof(long long uint) * (rep->maxLevel + 1);
        totalByte += rep->bl->getSize();
        totalByte += rep->bn->getSize();
        totalByte += rep->bt->getSize();
    }
    return totalByte;
}
