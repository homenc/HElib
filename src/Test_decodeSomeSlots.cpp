#include <NTL/ZZX.h>
#include "FHE.h"
#include "EncryptedArray.h"
#include "Ctxt.h"
#include "NumbTh.h"
int TestIt(long p, long m, long L)
{
    FHEcontext context(m, p, 1);
    buildModChain(context, L);
    FHESecKey sk(context);
    sk.GenSecKey(64);
    const FHEPubKey &pk = sk;
    auto ea = context.ea;
    std::vector<long> slots(ea->size());
    std::vector<NTL::ZZX> polys(ea->size());

    for (long i = 0; i < ea->size(); i++) {
        slots[i] = i + 1;
        NTL::SetCoeff(polys[i], 0, i + 1);
        NTL::SetCoeff(polys[i], 1, i + 1); // (i + 1) + (i + 1)X
    }

    Ctxt ctx1(pk), ctx2(pk);
    ea->encrypt(ctx1, pk, slots);
    ea->encrypt(ctx2, pk, polys);

    std::vector<long> positions;
    for (long i = 0; i < 4 && i < ea->size(); i++) {
        positions.push_back(i + 1);
    }

    NTL::ZZX ptxt;
    sk.Decrypt(ptxt, ctx1);
    ea->decodeSomeSlots(slots, ptxt, positions);
    for (size_t i = 0; i < positions.size(); i++) {
        assert(slots[i] == (positions[i] + 1) && "Fail at decodeSomeSlots<long>");
        printf("slot %ld is %ld\n", positions[i], slots[i]);
    }

    sk.Decrypt(ptxt, ctx2);
    ea->decodeSomeSlots(polys, ptxt, positions);
    for (size_t i = 0; i < positions.size(); i++) {
        assert(NTL::coeff(polys[i], 0) == positions[i] + 1 && "Fail at decodeSomeSlots<ZZX>");
        assert(NTL::coeff(polys[i], 1) == positions[i] + 1 && "Fail at decodeSomeSlots<ZZX>");
        std::cout << "slot " << positions[i] << " is " << polys[i] << "\n";
    }
    return 0;
}

int main(int argc, char *argv[]) {
    ArgMapping args;
    long p = 8191, m = 16384, L = 5;
    args.arg("p", p, "plaintext");
    args.arg("m", m, "m");
    args.arg("L", L, "level");
    args.parse(argc, argv);
    return TestIt(p, m, L);
}
