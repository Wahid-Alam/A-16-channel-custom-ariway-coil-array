function left = Eval_grad_left(f,A,At,opts)

f=reshape(f, opts.nx,opts.ny,opts.nz);
left = At(A(f));
left = left(:);

