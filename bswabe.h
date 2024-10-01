// bswabe.h
#ifndef BSWABE_H
#define BSWABE_H
#include <pbc/pbc.h>
#include <string>
#include <vector>

/***************************************part1: SetUp阶段生成系统公钥BswabePub和主私钥BswabeMsk***********************************************/
// 公钥结构体，包含系统配对和公钥元素
struct BswabePub {
    pairing_t pairing;   // 双线性配对参数
    element_t g;         // 群生成元 g ∈ G_0
    element_t h;         // h = g^β ∈ G_0
    element_t f;         // f = g^(1/β) ∈ G_0
    element_t e_gg_alpha;// e(g, g)^α ∈ G_T
};
// 主私钥结构体
struct BswabeMsk {
    element_t beta;      // β ∈ Z_p
    element_t g_alpha;   // g^α ∈ G_0
};

/***************************************part2: KeyGen阶段生成用户私钥：SK={D,D_j,D_j'}***********************************************/
// 用户私钥中的属性组件结构体
struct BswabePrvComp {
    std::string attr;    // 属性名
    element_t d;         // D_j = g^r * H(j)^r_j ∈ G_0
    element_t dp;        // D_j' = g^r_j ∈ G_0
};
// 用户私钥结构体
struct BswabePrv {
    element_t d;         // D = g^((α + r) / β)
    std::vector<BswabePrvComp*> comps; // 属性组件列表
};

/***************************************part3: Encrypt阶段生成密文：CT={BswabePolicy(树结构),cs,c,C_y,C_y'}***********************************************/
// 多项式结构体，用于加密策略树的每个节点
struct BswabePolynomial {
    int deg;                        // 多项式的度
    std::vector<element_s*> coef;   // 多项式系数
};
// 访问策略节点结构体
struct BswabePolicy {
    int k;                // 阈值 k
    std::string attr;     // 属性名（叶子节点）
    BswabePolynomial* q;  // 多项式 q_x
    std::vector<BswabePolicy*> children; // 子节点
    element_t c;          // 密文CT中C_y = g^{q_y(0)} ∈ G_0
    element_t cp;         // 密文CT中C_y' = H(att(y))^{q_y(0)} ∈ G_0
    bool satisfiable;     // 用于解密时是否满足条件
    std::vector<int> satl;// 满足条件的子节点索引列表: 从1开始...
};
// 密文结构体
struct BswabeCph {
    element_t cs;     // 加密后的明文 \tilde{C}
    element_t c;      // h^s
    BswabePolicy* p;  // 访问策略的根节点
};

/***************************************part4: Decrypt阶段生成明文：M={e,b}***********************************************/
// 解密结果结构体
struct BswabeElementBoolean {
    element_t e;  // 解密结果
    bool b;       // 解密成功标志
};


/***************************************函数声明***********************************************/
// 1.SetUp：系统初始化函数，生成公钥和主私钥
void bswabe_setup(BswabePub* pub, BswabeMsk* msk);

// 2.KeyGen：用户私钥生成函数，基于主私钥和用户属性生成用户私钥
BswabePrv* bswabe_keygen(BswabePub* pub, BswabeMsk* msk, const std::vector<std::string>& attrs);

// 3.Encrypt：加密函数，根据访问策略加密明文
BswabeCph* bswabe_enc(BswabePub* pub, const std::string& policy_str, element_t m);

// 4.Decrypt：解密函数，使用用户私钥解密密文
BswabeElementBoolean* bswabe_dec(BswabePub* pub, BswabePrv* prv, BswabeCph* cph);

// 释放用户私钥
void free_bswabe_prv(BswabePrv* prv);

// 释放密文
void free_bswabe_cph(BswabeCph* cph);

// 打印策略树（用于调试）
void print_policy_tree(BswabePolicy* p, int level = 0);

#endif // BSWABE_H
