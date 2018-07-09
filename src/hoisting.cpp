#include "hoisting.h"
#include "Ctxt.h"
#include "FHE.h"
#include <NTL/xdouble.h>

struct BasicAutomorphPrecon::Impl {
    Ctxt ctxt;
    NTL::xdouble noise;
    std::vector<DoubleCRT> polyDigits;

    Impl(const Ctxt& _ctxt) : ctxt(_ctxt), noise(1.0) {
        FHE_TIMER_START;
        if (ctxt.parts.size() >= 1) assert(ctxt.parts[0].skHandle.isOne());
        if (ctxt.parts.size() <= 1) return; // nothing to do

        ctxt.cleanUp();
        const FHEcontext& context = ctxt.getContext();
        const FHEPubKey& pubKey = ctxt.getPubKey();
        long keyID = ctxt.getKeyID();

        // The call to cleanUp() should ensure that this assertions passes.
        assert(ctxt.inCanonicalForm(keyID));

        // Compute the number of digits that we need and the esitmated
        // added noise from switching this ciphertext.
        long nDigits;
        std::tie(nDigits, noise)
            = ctxt.computeKSNoise(1, pubKey.keySWlist().at(0).ptxtSpace);

        double logProd = context.logOfProduct(context.specialPrimes);
        noise += ctxt.getNoiseVar() * xexp(2*logProd);

        // Break the ciphertext part into digits, if needed, and scale up these
        // digits using the special primes.

        ctxt.parts[1].breakIntoDigits(polyDigits, nDigits);
    }

    std::shared_ptr<Ctxt> automorph(long k) const {
        FHE_TIMER_START;

        // A hack: record this automorphism rather than actually performing it
        if (isSetAutomorphVals()) { // defined in NumbTh.h
            recordAutomorphVal(k);
            return make_shared<Ctxt>(ctxt);
        }

        if (k==1 || ctxt.isEmpty()) return make_shared<Ctxt>(ctxt);// nothing to do

        const FHEcontext& context = ctxt.getContext();
        const FHEPubKey& pubKey = ctxt.getPubKey();
        shared_ptr<Ctxt> result = make_shared<Ctxt>(ZeroCtxtLike, ctxt); // empty ctxt
        result->noiseVar = noise; // noise estimate

        if (ctxt.parts.size()==1) { // only constant part, no need to key-switch
            CtxtPart tmpPart = ctxt.parts[0];
            tmpPart.automorph(k);
            tmpPart.addPrimesAndScale(context.specialPrimes);
            result->addPart(tmpPart, /*matchPrimeSet=*/true);
            return result;
        }

        // Ensure that we have a key-switching matrices for this automorphism
        long keyID = ctxt.getKeyID();
        if (!pubKey.isReachable(k,keyID)) {
            throw std::logic_error("no key-switching matrices for k="+std::to_string(k)
                                   + ", keyID="+std::to_string(keyID));
        }

        // Get the first key-switching matrix for this automorphism
        const KeySwitch& W = pubKey.getNextKSWmatrix(k,keyID);
        long amt = W.fromKey.getPowerOfX();

        // Start by rotating the constant part, no need to key-switch it
        CtxtPart tmpPart = ctxt.parts[0];
        tmpPart.automorph(amt);
        tmpPart.addPrimesAndScale(context.specialPrimes);
        result->addPart(tmpPart, /*matchPrimeSet=*/true);

        // Then rotate the digits and key-switch them
        vector<DoubleCRT> tmpDigits = polyDigits;
        for (auto&& tmp: tmpDigits) // rotate each of the digits
            tmp.automorph(amt);

        result->keySwitchDigits(W, tmpDigits); // key-switch the digits

        long m = context.zMStar.getM();
        if ((amt-k)%m != 0) { // amt != k (mod m), more automorphisms to do
            k = MulMod(k, InvMod(amt,m), m); // k *= amt^{-1} mod m
            result->smartAutomorph(k);       // call usual smartAutomorph
        }
        return result;
    }

    std::shared_ptr<Ctxt> frobeniusAutomorph(long j) const {
        long m = ctxt.getContext().zMStar.getPhiM();
        long p = ctxt.getContext().zMStar.getP();
        long d = ctxt.getContext().zMStar.getOrdP();
        j = mcMod(j, d);
        return automorph(NTL::PowerMod(p, j, m));
    }
};


BasicAutomorphPrecon::BasicAutomorphPrecon(const Ctxt& _ctxt) 
    : impl_(std::make_shared<Impl>(_ctxt))
{ }

std::shared_ptr<Ctxt> BasicAutomorphPrecon::automorph(long k) const {
    if (impl_)
        return impl_->automorph(k);
    return nullptr;
}

std::shared_ptr<Ctxt> BasicAutomorphPrecon::frobeniusAutomorph(long j) const {
    if (impl_)
        return impl_->frobeniusAutomorph(j);
    return nullptr;
}

