#pragma once
#include <memory>
class Ctxt;
class BasicAutomorphPrecon {
public:
    BasicAutomorphPrecon(const Ctxt& _ctxt);

    void automorph(Ctxt &result, long k) const;
    void frobeniusAutomorph(Ctxt &result, long j) const;

    std::shared_ptr<Ctxt> automorph(long k) const;
    std::shared_ptr<Ctxt> frobeniusAutomorph(long j) const;
private:
    struct Impl;
    std::shared_ptr<Impl> impl_;
};
