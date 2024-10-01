#include "bswabe.h"
#include <openssl/sha.h>
#include <iostream>
#include <sstream>
#include <stack>
#include <map>
#include <cstring>
#include <algorithm>    // 为 std::find 提供支持
#include <sstream>      // 引入字符串流库

// 哈希函数，将字符串映射为群元素
void element_from_string(element_t h, const std::string& s) {
    unsigned char digest[SHA_DIGEST_LENGTH];
    SHA1(reinterpret_cast<const unsigned char*>(s.c_str()), s.length(), digest);
    element_from_hash(h, digest, SHA_DIGEST_LENGTH);  // 将哈希值转换为群元素
}

/**
 * 1.SetUp:系统初始化函数
 *
 * 功能描述：
 * 此函数用于初始化基于双线性对的加密系统，生成公共参数（公钥）和主私钥。它读取配对参数文件，初始化配对环境，并生成系统所需的公共元素和私有元素。
 *
 * 输入参数：
 * - pub：BswabePub*，公钥结构体指针，将在函数中被初始化和设置。
 * - msk：BswabeMsk*，主私钥结构体指针，将在函数中被初始化和设置。
 *
 * 返回值：
 * - 无。函数通过修改传入的 pub 和 msk 指针来设置公钥和主私钥。
 *
 * 注意事项：
 * - 在调用此函数之前，pub 和 msk 需要已分配内存，但其中的元素无需初始化。
 * - 函数内部会读取配对参数文件，路径需要根据实际环境修改，确保参数文件存在并可读取。
 * - 函数结束后，pub 和 msk 中的元素已初始化，使用完毕后需要调用适当的清理函数释放内存。
 * - 函数内部使用的临时元素在函数结束前已被正确清理。
 * - 如果无法打开参数文件或读取失败，函数将打印错误信息并退出程序。
 */
void bswabe_setup(BswabePub* pub, BswabeMsk* msk) {
    element_t alpha, beta_inv;
    pbc_param_t param;
    // 初始化pairing，通过从参数文件读取
    char param_buf[2048];
    FILE* param_file = fopen("/home/hututu/instllpkg/pbc-0.5.14/param/a.param", "r"); // 替换为实际的参数文件路径
    if (!param_file) {
        fprintf(stderr, "Error opening param file\n");
        exit(1);
    }
    // 读取文件内容
    size_t count = fread(param_buf, 1, 2048, param_file);
    fclose(param_file);
    if (!count) pbc_die("Error reading param file");
    // 根据读取的参数初始化 pairing
    pairing_init_set_buf(pub->pairing, param_buf, count);
    // 初始化公钥和私钥的元素
    element_init_G1(pub->g, pub->pairing);
    element_init_G1(pub->h, pub->pairing);
    element_init_G1(pub->f, pub->pairing);
    element_init_GT(pub->e_gg_alpha, pub->pairing);
    element_init_Zr(alpha, pub->pairing);
    element_init_Zr(msk->beta, pub->pairing);
    element_init_G1(msk->g_alpha, pub->pairing);
    // 随机选择 α 和 β
    element_random(alpha);      // 随机选择 α ∈ Z_p
    element_random(msk->beta);  // 随机选择 β ∈ Z_p
    element_random(pub->g);     // 随机生成元 g ∈ G1
    // 计算 g^α
    element_pow_zn(msk->g_alpha, pub->g, alpha);  // g^α
    // 计算 β 的逆元并使用它进行相关计算
    element_init_Zr(beta_inv, pub->pairing);      // β 的逆元
    element_invert(beta_inv, msk->beta);          // beta_inv = 1/β
    element_pow_zn(pub->f, pub->g, beta_inv);     // f = g^{1/β}
    element_pow_zn(pub->h, pub->g, msk->beta);    // h = g^β
    // 计算 e(g, g)^α
    element_pairing(pub->e_gg_alpha, pub->g, msk->g_alpha); // e(g, g)^α
    // 清除临时元素
    element_clear(alpha);
    element_clear(beta_inv);
}

/**
 * 2.KeyGen:用户私钥生成函数
 *
 * 功能描述：
 * 根据系统公钥、公钥参数和用户的属性集合，生成用户的私钥。私钥包括主私钥部分 D，以及与用户属性相关的部分 D_j 和 D_j'。
 *
 * 输入参数：
 * - pub：BswabePub*，系统公钥结构体，包含公共参数（如 g、h、f、e_gg_alpha 等）。
 * - msk：BswabeMsk*，主私钥结构体，包含私有参数（如 β、g^α 等）。
 * - attrs：const std::vector<std::string>&，用户属性的字符串列表。
 *
 * 返回值：
 * - BswabePrv*，生成的用户私钥结构体，包含 D 以及属性组件 D_j 和 D_j'。
 *
 * 注意事项：
 * - 确保在调用此函数之前，公钥 pub 和主私钥 msk 已正确初始化。
 * - 返回的 BswabePrv* 指针需要在使用完毕后手动释放，防止内存泄漏。
 * - 函数内部使用的临时元素在函数结束前已被正确清理。
 * - 此函数不验证输入参数的合法性，调用者需确保输入的正确性。
 */
BswabePrv* bswabe_keygen(BswabePub* pub, BswabeMsk* msk, const std::vector<std::string>& attrs) {
    BswabePrv* prv = new BswabePrv();    // 用户私钥，D,D_j,D_j'
    element_t r,g_r,beta_inv;            // 随机数r ∈ Z_p
    
    element_init_Zr(r, pub->pairing);
    element_init_G1(g_r, pub->pairing);
    element_init_G1(prv->d, pub->pairing); 
    element_init_Zr(beta_inv, pub->pairing);
    // 1. 随机选择 r ∈ Z_p
    element_random(r);  // 随机选择 r ∈ Z_p
    // 2. 计算 beta 的逆元 beta_inv = 1 / beta
    element_invert(beta_inv, msk->beta);
    // 3. 计算 g^r
    element_pow_zn(g_r, pub->g, r);    // g_r = g^r
    // 4. 计算 g^{alpha + r} = g^alpha * g^r
    element_mul(prv->d, msk->g_alpha, g_r);     // D = g^alpha * g^r = g^(α + r)
    // 5. 计算 D = (g^{alpha + r})^{beta_inv}, 现在 D 即为所需的私钥元素 temp = g^{(α + r)/β}
    element_pow_zn(prv->d, prv->d, beta_inv);    // D = temp^{beta_inv}
    // 6. 生成与用户属性相关的：D_j, D_j'
    for (const std::string& attr : attrs) {
        BswabePrvComp* comp = new BswabePrvComp();
        comp->attr = attr;

        element_t rj;
        element_init_Zr(rj, pub->pairing);
        element_random(rj);                         // 随机选择 r_j ∈ Z_p
        element_init_G1(comp->d, pub->pairing);     // D_j
        element_init_G1(comp->dp, pub->pairing);    // D_j'
        element_t h_attr,h_attr_rj;
        element_init_G1(h_attr, pub->pairing);
        element_init_G1(h_attr_rj, pub->pairing);
        element_from_string(h_attr, comp->attr);    // 根据用户属性string生成H(j)
        element_pow_zn(h_attr_rj, h_attr, rj);      // h_attr_rj = H(j)^rj
        element_pow_zn(comp->d, pub->g, r);         // D_j = g^r
        element_mul(comp->d, comp->d, h_attr_rj);   // g^r * H(j)^r_j
        element_pow_zn(comp->dp, pub->g, rj);       // D_j' = g^{r_j}
        prv->comps.push_back(comp);                 // 将与某一用户属性相关的D_j与D_j'存入prv->comps中
        element_clear(h_attr);
        element_clear(h_attr_rj);
        element_clear(rj);
    }

    // 清除元素
    element_clear(r);
    element_clear(beta_inv);
    element_clear(g_r);
    return prv;
}

// 生成访问策略的节点
BswabePolicy* base_node(int k, const std::string& s) {
    BswabePolicy* p = new BswabePolicy();
    p->k = k;
    p->attr = s;
    p->q = nullptr;
    p->satisfiable = false;
    return p;
}

/**
 * 解析访问策略的后缀表达式并生成策略树函数
 *
 * 功能描述：此函数用于解析给定的后缀表达式形式的访问策略字符串，并将其转换为策略树（`BswabePolicy` 结构体）。策略字符串采用逆波兰表示法，其中属性和阈值门以空格分隔。函数通过遍历令牌并使用栈来构建策略树。
 *
 * 输入参数：
 * - `s`：`const std::string&`，后缀表达式形式的策略字符串，例如 `"attr1 attr2 1of2"`。
 *
 * 返回值：
 * - `BswabePolicy*`，解析得到的策略树的根节点指针。如果解析失败，程序将输出错误信息并终止。
 *
 * 注意事项：
 * - 策略字符串应按正确的后缀表达式格式编写，属性和阈值门之间以空格分隔。
 * - 阈值门的格式为 `"kofn"`，表示在 `n` 个子节点中至少满足 `k` 个条件。
 * - 函数假定输入的策略字符串格式正确，不进行格式验证。
 * - 函数使用了 `exit(1)` 在发生错误时终止程序，实际使用中可修改为抛出异常或其他错误处理方式。
 */
BswabePolicy* parse_policy_postfix(const std::string& s) {
    std::vector<std::string> tokens;            // 用于存储拆分后的策略令牌
    std::istringstream iss(s);
    std::string token;
    // 将策略字符串按空格拆分为令牌，并存入 tokens 向量
    while (iss >> token) {
        tokens.push_back(token);
    }
    std::vector<BswabePolicy*> stack;           // 栈，用于构建策略树

    // 遍历每个令牌，构建策略树
    for (const auto& tok : tokens) {
        if (tok.find("of") == std::string::npos) {
            // 如果令牌不包含 "of"，则认为是属性，创建叶子节点并压入栈
            stack.push_back(base_node(1, tok));     // 叶子节点，阈值 k=1，属性名为 tok
        } else {
            // 令牌包含 "of"，表示是阈值门，需要解析 k 和 n
            size_t pos = tok.find("of");
            int k = std::stoi(tok.substr(0, pos));  // 提取阈值 k
            int n = std::stoi(tok.substr(pos + 2)); // 提取子节点数量 n
            // 检查栈中是否有足够的节点供弹出
            if (stack.size() < n) {
                std::cerr << "错误：栈中节点不足，无法弹出所需的子节点数" << std::endl;
                exit(1);
            }   
            // 创建新的策略节点，阈值为 k，属性名为空（因为是内部节点）
            BswabePolicy* node = base_node(k, "");  // 阈值节点，阈值 k=k，属性名为空
            node->children.resize(n);               // 调整子节点向量的大小
            // 从栈中弹出 n 个节点，作为当前节点的子节点
            for (int i = n - 1; i >= 0; --i) {
                node->children[i] = stack.back();   // 从栈顶获取子节点
                stack.pop_back();                   // 弹出栈顶元素
            }
            // 将新创建的节点压入栈
            stack.push_back(node);
        }
    }
    // 解析结束后，栈中应只剩下一个节点，即根节点
    if (stack.size() != 1) {
        std::cerr << "错误：策略字符串格式不正确，无法生成唯一的根节点" << std::endl;
        exit(1);
    }
    return stack.back(); // 返回根节点
}

// 打印策略树（用于调试）
void print_policy_tree(BswabePolicy* p, int level) {
    for (int i = 0; i < level; ++i)
        std::cout << "  ";
    std::cout << "Node: k = " << p->k << ", attr = " << p->attr << std::endl;
    for (const auto& child : p->children) {
        print_policy_tree(child, level + 1);
    }
}

/**
 * 生成随机多项式函数
 *
 * 功能描述：生成一个指定度数的随机多项式，并使该多项式在零点处的值为指定的 `zero_val`。多项式用于加密策略树中的各个节点，以实现访问控制策略。
 *
 * 输入参数：
 * - deg：int，多项式的度数。
 * - zero_val：element_t，指定的零点处的值，即多项式在 x=0 处的值。
 *
 * 返回值：
 * - BswabePolynomial*，生成的多项式结构体指针，包含多项式的度数和系数。
 *
 * 注意事项：
 * - 返回的多项式结构体需要在使用完毕后手动释放其内存，包括多项式系数中的元素，防止内存泄漏。
 * - 函数内部使用了 PBC 库的元素类型，调用者需确保 PBC 环境已正确初始化。
 * - 为了防止内存泄漏，函数在每次分配元素后，都需要确保在程序结束时对其进行释放。
 */
BswabePolynomial* rand_poly(int deg, element_t zero_val) {
    BswabePolynomial* q = new BswabePolynomial(); // 创建新的多项式结构体
    q->deg = deg;                                 // 设置多项式的度数
    q->coef.resize(deg + 1);                      // 调整系数向量的大小，系数数量为 deg + 1
    // 初始化多项式的系数元素
    for (int i = 0; i <= deg; i++) {
        q->coef[i] = (element_s*)malloc(sizeof(element_s));     // 给每个系数coef分配存储空间
        element_init_same_as(q->coef[i], zero_val);             // 初始化每个系数元素，与 zero_val 类型相同
    }
    // 设置多项式在 x=0 处的值为 zero_val，即常数项，也即s=qx(0)
    element_set(q->coef[0], zero_val);
    // 为多项式的其他系数赋值
    for (int i = 1; i <= deg; i++) {
        element_random(q->coef[i]);                 // 随机生成系数元素
        // element_set_si(q->coef[i], 2);           // （可选）将系数设置为固定值，例如 2，用于测试
    }
    return q; // 返回生成的多项式结构体指针
}

/**
 * 评估多项式在指定点的值函数
 *
 * 功能描述：计算给定多项式 `q` 在点 `x` 处的值，并将结果存储在元素 `r` 中。多项式使用系数表示，其中系数为群元素类型（`element_t`）。该函数适用于基于双线性对的密码系统中的多项式求值操作。
 *
 * 输入参数：
 * - `q`：`BswabePolynomial*`，指向多项式结构体的指针，包含多项式的度数和系数。
 * - `x`：`element_t`，多项式求值的点。调用前需要初始化。
 * 
 * 输出参数：
 * - `r`：`element_t`，用于存储计算结果的元素。调用前需要初始化。
 * 
 * 返回值：
 * - 无。计算结果直接存储在参数 `r` 中。
 *
 * 注意事项：
 * - 在调用此函数之前，参数 `r`、`x` 和多项式 `q` 的系数都必须已正确初始化，且与配对环境匹配。
 * - 函数内部使用了临时元素 `sum`、`exp`、`term`，在函数结束前已被正确清理，防止内存泄漏。
 * - 调用者需确保多项式 `q` 的度数和系数已正确设置。
 * - 该函数假定多项式的系数存储在 `q->coef` 中，索引从 0 到 `q->deg`。
 */
void eval_poly(element_t r, BswabePolynomial* q, element_t x) {
    element_t sum, exp, term;
    // 初始化临时变量 sum，用于累加多项式求和值，与 r 类型相同
    element_init_same_as(sum, r);
    // 初始化临时变量 term，用于存储每一项的计算结果，与 r 类型相同
    element_init_same_as(term, r);
    // 初始化临时变量 exp，用于计算 x 的幂次，与 x 类型相同
    element_init_same_as(exp, x);
    // sum = 0，初始化求和值为零
    element_set0(sum);
    // exp = 1，初始化 x 的指数次幂为 1（即 x^0）
    element_set1(exp);
    // 遍历多项式的所有系数，计算 sum = ∑ (q->coef[i] * x^i)
    // 例如：f(x) = ax^2+bx+c，一个三个系数，下面循环执行三次；第一次结束sum = c, exp = x; 第二次结束sum = c + bx, exp = x^2; 第三次结束sum = c + bx + ax^2, exp = x^3;
    for (int i = 0; i <= q->deg; ++i) {
        // term = q->coef[i] * exp，计算当前项的值
        element_mul(term, q->coef[i], exp);
        // sum += term，将当前项的值累加到总和中
        element_add(sum, sum, term);
        // exp *= x，更新 exp 为 x 的下一个幂次（即 x^{i+1}）
        element_mul(exp, exp, x);
    }
    // 将计算得到的多项式值 sum 赋值给结果变量 r，返回的结果变量r
    element_set(r, sum);
    // 清理临时变量，释放内存
    element_clear(sum);
    element_clear(exp);
    element_clear(term);
}

/**
 * 填充策略树并计算加密组件函数，主要用在加密过程，对叶子节点（整个树的所有叶子节点）进行加密
 *
 * 功能描述：该函数用于在加密过程中，递归地填充策略树的每个节点。对于每个节点，生成一个随机多项式 `q`，并根据该多项式计算节点的加密组件。对于叶子节点，计算并存储加密所需的群元素；对于内部节点，递归处理其子节点。
 *
 * 输入参数：
 * - `p`：`BswabePolicy*`，指向策略树节点的指针，需要填充的节点。
 * - `pub`：`BswabePub*`，系统公钥结构体，包含公共参数（如 `g`、`h`、`pairing` 等）。
 * - `e`：`element_t`，上级节点在当前节点处的多项式值，即 `q(0)` 的值。对于根节点，`e` 通常为随机生成的秘密值 `s`。
 *
 * 返回值：
 * - 无。函数通过修改参数 `p` 来填充策略树节点的内容。
 *
 * 注意事项：
 * - 调用此函数前，策略树应已正确构建（例如通过 `parse_policy_postfix` 函数）。
 * - 函数内部使用了递归调用来处理策略树的所有节点。
 * - 函数中使用了临时元素，需要确保在函数结束前正确清理，防止内存泄漏。
 * - 对于叶子节点，函数会初始化并计算加密所需的元素 `p->c` 和 `p->cp`。
 */
void fill_policy(BswabePolicy* p, BswabePub* pub, element_t e) {
    // 生成随机多项式 q，并设置其在 x=0 处的值为 e
    p->q = rand_poly(p->k - 1, e);  // 多项式的度为 k-1，q(0) = e

    if (p->children.empty()) {      // 处理叶子节点
        // 初始化叶子节点的加密组件 c 和 cp
        element_init_G1(p->c, pub->pairing);
        element_init_G1(p->cp, pub->pairing);
        // 将属性字符串映射为群元素 h_attr = H(attr)
        element_t h_attr;
        element_init_G1(h_attr, pub->pairing);
        element_from_string(h_attr, p->attr);
        // 计算加密组件
        // c = g^{q(0)}
        element_pow_zn(p->c, pub->g, p->q->coef[0]);
        // cp = H(attr)^{q(0)}
        element_pow_zn(p->cp, h_attr, p->q->coef[0]);
        // 清理临时元素 h_attr
        element_clear(h_attr);
    } else {        // 处理阈值节点
        // 初始化临时元素 index，用于表示子节点的序号
        element_t index;
        element_init_Zr(index, pub->pairing);
        // 遍历子节点
        for (size_t i = 0; i < p->children.size(); ++i) {
            // 设置 index = i + 1，因为子节点序号从 1 开始
            element_set_si(index, i + 1);
            // 初始化临时元素 q_y0，用于存储多项式在子节点序号处的值 q(i+1)
            element_t q_y0;
            element_init_Zr(q_y0, pub->pairing);
            // 计算多项式 q 在 index 处的值，即 f(x) = q_y0(父节点的f(x)值为子节点的秘密值) = q(index)
            eval_poly(q_y0, p->q, index);
            // 递归调用 fill_policy，填充子节点，传递 q_y0 作为新的 e 值
            fill_policy(p->children[i], pub, q_y0);
            // 清理临时元素 q_y0
            element_clear(q_y0);
        }
        // 清理临时元素 index
        element_clear(index);
    }
}

/**
 * 3.Encrypt:加密函数
 *
 * 功能描述：根据给定的访问策略字符串，对明文元素进行加密。函数生成一个新的密文结构体 `BswabeCph`，其中包含加密的明文和访问策略树。加密过程中使用随机生成的秘密值 `s`，并将其与访问策略树结合，生成策略相关的加密组件。
 *
 * 输入参数：
 * - `pub`：`BswabePub*`，系统公钥结构体，包含公共参数（如 `g`、`h`、`e_gg_alpha` 等）。
 * - `policy_str`：`const std::string&`，访问策略字符串，采用后缀表达式形式，例如 `"attr1 attr2 1of2"`。
 * - `m`：`element_t`，要加密的明文元素，属于 `GT` 群。
 *
 * 返回值：
 * - `BswabeCph*`，生成的密文结构体，包含加密的明文和策略树。
 *
 * 注意事项：
 * - 在调用此函数之前，公钥 `pub` 必须已正确初始化。
 * - 返回的密文结构体需要在使用完毕后手动释放，防止内存泄漏。
 * - 函数内部生成的随机秘密值 `s` 是加密过程中的关键参数，确保其随机性和安全性。
 * - 函数中使用的临时元素在函数结束前已被正确清理。
 */
BswabeCph* bswabe_enc(BswabePub* pub, const std::string& policy_str, element_t m) {
    BswabeCph* cph = new BswabeCph();           // 创建新的密文结构体
    element_t s;                                // 随机秘密值 s，用于加密过程
    element_init_Zr(s, pub->pairing);           // 初始化元素 s，属于 Z_p 群
    element_random(s);                          // 生成随机秘密值 s ∈ Z_p  
    // 初始化密文结构体中的元素
    element_init_GT(cph->cs, pub->pairing);    // 初始化加密后的明文元素 cs，属于 GT 群
    element_init_G1(cph->c, pub->pairing);     // 初始化加密组件 c，属于 G1 群

    /***********************************CT第一部分：C_y,C_y'************************************************/
    // 解析访问策略字符串，生成策略树
    cph->p = parse_policy_postfix(policy_str); // 解析访问策略，生成策略树 cph->p
    // 填充策略树，生成策略相关的加密组件
    fill_policy(cph->p, pub, s);               // 填充策略树，使用秘密值 s

    /***********************************CT第二部分：C=h^s, C_head=Me(g,g)^{α*s}************************************************/
    // 计算加密的明文部分
    element_t egg_alpha_s;                          // 临时元素，用于存储 e(g, g)^{αs}
    element_init_GT(egg_alpha_s, pub->pairing);
    // 计算 egg_alpha_s = e(g, g)^{αs}
    element_pow_zn(egg_alpha_s, pub->e_gg_alpha, s); // egg_alpha_s = (e(g, g)^α)^s = e(g, g)^{αs}
    // 计算密文中的加密明文部分 cph->cs = m * e(g, g)^{αs}
    element_mul(cph->cs, m, egg_alpha_s);            // C_head = cs = m * e(g, g)^{αs}
    // 计算密文中的加密组件 c = h^s
    element_pow_zn(cph->c, pub->h, s);               // c = h^s

    // 清理临时元素
    element_clear(s);                                // 清理秘密值 s
    element_clear(egg_alpha_s);                      // 清理临时元素 egg_alpha_s

    return cph; // 返回生成的密文结构体
}


/**
 * 检查策略是否满足并标记满足条件的节点函数
 *
 * 功能描述：此函数递归地检查策略树中的每个节点，确定用户的属性集合是否满足策略要求。
 * 对于叶子节点，函数检查其属性是否在用户的属性集合中。
 * 对于中间节点，函数统计满足条件的子节点数量，并根据阈值 `k` 决定当前节点是否满足。
 * 函数会设置节点的 `satisfiable` 标志，并在中间节点的 `satl` 列表中记录满足条件的子节点索引。
 * 该函数在解密过程中用于确定用户是否有足够的属性来满足加密策略，从而能够成功解密密文。
 *
 * 输入参数：
 * - `p`：`BswabePolicy*`，指向策略树节点的指针，需要检查的节点。
 * - `attrs`：`const std::vector<std::string>&`，用户持有的属性集合。
 *
 * 返回值：
 * - `bool`，表示当前节点是否满足策略要求。`true` 表示满足，`false` 表示不满足。
 *
 * 注意事项：
 * - 函数通过递归调用，遍历整个策略树。
 * - 对于叶子节点，函数会检查其属性是否在用户属性集合中，并设置 `satisfiable` 标志。
 * - 对于中间节点，函数会递归检查其子节点，统计满足条件的子节点数量，并根据阈值 `k` 决定当前节点是否满足。
 * - 函数会修改策略树节点的状态，包括 `satisfiable` 标志和 `satl` 列表。
 * - 在解密过程中，需要在调用此函数后，才能正确地进行后续的拉格朗日插值等操作。
 */
bool check_sat(BswabePolicy* p, const std::vector<std::string>& attrs) {
    if (p->children.empty()) { // 叶子节点处理
        // 检查叶子节点的属性是否在用户的属性集合中
        // 如果存在，则设置 satisfiable 为 true；否则为 false
        p->satisfiable = std::find(attrs.begin(), attrs.end(), p->attr) != attrs.end();
    } else {        // 阈值节点处理
        int satisfied = 0;           // 统计满足条件的子节点数量
        p->satl.clear();             // 清空满足条件的子节点索引列表
        // 遍历所有子节点，递归检查每个子节点是否满足条件
        for (size_t i = 0; i < p->children.size(); ++i) {
            if (check_sat(p->children[i], attrs)) {
                satisfied++;                // 如果子节点满足条件，增加计数
                p->satl.push_back(i + 1);   // 记录满足条件的子节点索引，索引从 1 开始计数
            }
        }
        // 根据阈值 k，判断当前节点是否满足条件
        // 如果满足条件的子节点数量不少于阈值 k，则当前节点满足
        p->satisfiable = (satisfied >= p->k);
    }
    // 返回当前节点的 satisfiable 状态
    return p->satisfiable;
}

/**
 * 拉格朗日插值系数计算函数
 *
 * 功能描述：该函数用于计算在拉格朗日插值中用于秘密共享重构的拉格朗日系数 λ_i。
 * 给定满足条件的节点索引列表 `satl`，以及特定的节点索引 `i`，
 * 计算拉格朗日系数 λ_i，其公式为：
 *
 *     λ_i = ∏_{j ∈ satl, j ≠ i} (-j) / (i - j)
 *
 * 该系数用于在秘密共享重构中，根据满足条件的子节点计算父节点的值。
 *
 * 输入参数：
 * - `satl`：`const std::vector<int>&`，满足条件的子节点索引列表。
 * - `i`：`int`，当前计算的节点索引。
 * - `pub`：`BswabePub*`，系统公钥结构体，包含配对参数。
 *
 * 输出参数：
 * - `coef`：`element_t`，用于存储计算结果的元素。调用前需要初始化。
 * 
 * 返回值：
 * - 无。计算结果直接存储在参数 `coef` 中。
 *
 * 注意事项：
 * - 参数 `coef` 在调用前必须已使用 `element_init_Zr` 初始化。
 * - 函数内部使用了临时元素 `num`、`denom` 和 `tmp`，在函数结束前已被正确清理。
 * - 索引 `i` 和 `j` 应该从 1 开始，与策略树中子节点的编号一致。
 * - 该函数假定 `satl` 中的索引都是正整数，且包含 `i`。
 */
void lagrange_coefficient(element_t coef, const std::vector<int>& satl, int i, BswabePub* pub) {
    element_t num, denom; // 分别用于存储分子和分母的中间计算结果
    element_init_Zr(num, pub->pairing);   // 初始化分子 num 为整数域元素
    element_init_Zr(denom, pub->pairing); // 初始化分母 denom 为整数域元素
    element_set1(num);    // num = 1，初始化分子为 1
    element_set1(denom);  // denom = 1，初始化分母为 1

    // 计算拉格朗日系数 λ_i = ∏_{j ≠ i} (-j) / (i - j)
    for (int j : satl) {
        if (j == i) {
            continue;  // 跳过当前节点索引 i
        }
        element_t tmp; // 临时变量用于存储中间结果
        element_init_Zr(tmp, pub->pairing);

        // 计算分子部分 num *= -j
        element_set_si(tmp, -j);          // tmp = -j
        element_mul(num, num, tmp);       // num = num * tmp

        // 计算分母部分 denom *= (i - j)
        element_set_si(tmp, i - j);       // tmp = i - j
        element_mul(denom, denom, tmp);   // denom = denom * tmp

        element_clear(tmp); // 清理临时变量 tmp
    }

    // 计算 denom 的逆元
    element_invert(denom, denom);  // denom = 1 / denom
    // 计算拉格朗日系数 coef = num * denom
    element_mul(coef, num, denom); // coef = num * denom

    // 清理分子和分母的中间结果
    element_clear(num);
    element_clear(denom);
}

/**
 * 使用拉格朗日插值法解密节点函数
 *
 * 功能描述：此函数递归地解密策略树中的每个节点，利用拉格朗日插值法在满足策略的情况下重构秘密。
 * 对于叶子节点，函数计算解密元素；对于内部节点，函数递归地处理满足条件的子节点，计算拉格朗日系数，重构父节点的值。
 * 函数最终在根节点重构出加密过程中使用的秘密值，用于解密密文。
 *
 * 输入参数：
 * - `p`：`BswabePolicy*`，指向策略树节点的指针，需要解密的节点。
 * - `prv`：`BswabePrv*`，用户的私钥结构体，包含用户的属性和对应的私钥组件。
 * - `cph`：`BswabeCph*`，密文结构体，包含加密的明文和策略树。
 * - `pub`：`BswabePub*`，系统公钥结构体，包含公共参数和配对。
 *
 * 输出参数：
 * - `r`：`element_t`，用于存储计算结果的元素。调用前需要初始化。
 * 
 * 返回值：
 * - 无。计算结果直接存储在参数 `r` 中。
 *
 * 注意事项：
 * - 在调用此函数之前，需要确保策略树已使用 `check_sat` 函数标记了满足条件的节点。
 * - 函数通过递归调用，遍历满足条件的策略树节点，计算用于解密的元素。
 * - 函数内部使用了临时元素，需要在函数结束前正确清理，防止内存泄漏。
 * - 索引 `i` 在 `p->satl` 中是从 1 开始的，因此在访问子节点时需要减 1。
 */
void decrypt_node_with_lagrange(element_t r, BswabePolicy* p, BswabePrv* prv, BswabeCph* cph, BswabePub* pub) {
    if (p->children.empty()) {        // 叶子节点的解密操作
        // 在用户私钥中查找与当前叶子节点属性匹配的私钥组件
        for (auto& comp : prv->comps) {
            if (comp->attr == p->attr) {  // 找到用户私钥对应的属性密文
                // 初始化临时元素 e1 和 e2，用于存储配对运算结果
                element_t e1, e2;
                element_init_GT(e1, pub->pairing);
                element_init_GT(e2, pub->pairing);
                // 计算 e1 = e(C_y, D_j)
                element_pairing(e1, p->c, comp->d);     // e(C_y, D_j)
                // 计算 e2 = e(C_y', D_j')
                element_pairing(e2, p->cp, comp->dp);   // e(C_y', D_j')
                // 计算 e2 的逆元 e2^{-1}
                element_invert(e2, e2);                 // e2^{-1}
                // 计算 r = e1 * e2 = e(C_y, D_j) / e(C_y', D_j')
                element_mul(r, e1, e2);                 // r = e1 * e2
                // 清理临时元素 e1 和 e2
                element_clear(e1);
                element_clear(e2);
                break; // 找到匹配的属性，退出循环（一次仅处理一个叶子节点）
            }
        }
    } else {        // 处理非叶子节点，使用拉格朗日插值还原秘密 s
        element_t Fx, t;
        // 初始化 Fx，用于累积子节点的解密结果
        element_init_GT(Fx, pub->pairing);
        element_set1(Fx);  // 初始化 Fx 为 1
        // 初始化 t，用于存储拉格朗日系数
        element_init_Zr(t, pub->pairing);
        // 遍历满足条件的子节点索引列表 p->satl
        for (int i : p->satl) {
            element_t share;
            // 初始化 share，用于存储子节点的解密结果
            element_init_GT(share, pub->pairing);
            // 递归解密子节点，注意子节点数组从 0 开始，而 satl 索引从 1 开始
            decrypt_node_with_lagrange(share, p->children[i - 1], prv, cph, pub);
            // 计算拉格朗日系数 t = λ_i
            lagrange_coefficient(t, p->satl, i, pub);
            // 计算 share = share^{λ_i}
            element_pow_zn(share, share, t);
            // 累积计算 Fx = Fx * share
            element_mul(Fx, Fx, share);
            // 清理临时元素 share
            element_clear(share);
        }
        // 将累积结果 Fx 赋值给输出参数 r
        element_set(r, Fx);
        // 清理临时元素 Fx 和 t
        element_clear(Fx);
        element_clear(t);
    }
}
/**
 * 4.Decrypt:解密函数
 *
 * 功能描述：该函数用于对属性基加密的密文进行解密。它首先检查用户的属性集合是否满足密文的访问策略，
 * 如果满足，则通过递归方式计算必要的中间值，最终恢复出原始的明文元素。
 *
 * 输入参数：
 * - `pub`：`BswabePub*`，系统公钥结构体，包含公共参数和配对环境。
 * - `prv`：`BswabePrv*`，用户的私钥结构体，包含用户的属性和对应的私钥组件。
 * - `cph`：`BswabeCph*`，密文结构体，包含加密的明文和访问策略树。
 *
 * 返回值：
 * - `BswabeElementBoolean*`，包含解密得到的明文元素和解密成功与否的标志。
 *   - `result->e`：`element_t`，解密得到的明文元素，属于目标群 `GT`。
 *   - `result->b`：`bool`，表示解密是否成功。`true` 表示成功，`false` 表示失败。
 *
 * 注意事项：
 * - 在调用此函数之前，需要确保公钥、私钥和密文都已正确初始化。
 * - 函数内部会检查用户的属性集合是否满足密文的访问策略，如果不满足，将返回解密失败的结果。
 * - 函数内部使用了临时元素，需要在函数结束前正确清理，防止内存泄漏。
 * - 返回的结果结构体需要在使用完毕后手动释放，并清理内部的元素。
 */
BswabeElementBoolean* bswabe_dec(BswabePub* pub, BswabePrv* prv, BswabeCph* cph) {
    BswabeElementBoolean* result = new BswabeElementBoolean();
    // 初始化结果元素 result->e，属于目标群 GT
    element_init_GT(result->e, pub->pairing);
    // 提取用户私钥中的属性列表，存入 attrs
    std::vector<std::string> attrs;    // 从私钥中提取用户属性集
    for (size_t i = 0; i < prv->comps.size(); ++i) {
        attrs.push_back(prv->comps[i]->attr);  // 提取每个组件的属性
    }
    // 检查用户的属性是否满足密文的访问策略
    if (!check_sat(cph->p, attrs)) {
        std::cout << "用户属性不满足策略要求" << std::endl;
        result->b = false; // 解密失败
        return result;
    }
    // 执行解密过程，首先通过递归解密得到 e(g, g)^{rs}
    element_t e_gg_rs;
    element_init_GT(e_gg_rs, pub->pairing);
    element_set1(e_gg_rs); // 初始化 e_gg_rs 为 1
    // 递归解密策略树，计算 e_gg_rs = e(g, g)^{rs}
    decrypt_node_with_lagrange(e_gg_rs, cph->p, prv, cph, pub);

    // 计算 e(C, D)
    element_t e_CD;
    element_init_GT(e_CD, pub->pairing);
    element_pairing(e_CD, cph->c, prv->d); // e_CD = e(C, D)
    // 计算 e(C, D) / e(g, g)^{rs}
    element_div(e_CD, e_CD, e_gg_rs); // e_CD = e(C, D) / e(g, g)^{rs}
    // 计算最终的明文 M = \tilde{C} / (e(C, D) / e(g, g)^{rs})
    element_div(result->e, cph->cs, e_CD); // result->e = \tilde{C} / e_CD
    // 清理临时变量
    element_clear(e_gg_rs);
    element_clear(e_CD);
    result->b = true; // 解密成功
    return result;
}

// 释放策略树资源
void free_bswabe_policy(BswabePolicy* p) {
    if (!p) return;
    for (auto& child : p->children) {
        free_bswabe_policy(child);
    }
    delete p;
}
// 释放用户私钥资源
void free_bswabe_prv(BswabePrv* prv) {
    if (!prv) return;
    for (auto& comp : prv->comps) {
        element_clear(comp->d);
        element_clear(comp->dp);
        delete comp;
    }
    delete prv;
}
// 释放密文资源
void free_bswabe_cph(BswabeCph* cph) {
    if (!cph) return;
    free_bswabe_policy(cph->p);
    element_clear(cph->cs);
    element_clear(cph->c);
    delete cph;
}
