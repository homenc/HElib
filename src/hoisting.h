#pragma once
#include <memory>
class Ctxt;
class BasicAutomorphPrecon {
public:
    BasicAutomorphPrecon(const Ctxt& _ctxt);

    std::shared_ptr<Ctxt> automorph(long k) const;
    std::shared_ptr<Ctxt> frobeniusAutomorph(long j) const;
private:
    struct Impl;
    std::shared_ptr<Impl> impl_;
};
