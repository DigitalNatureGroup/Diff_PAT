module Adamopt

# This is a module implementing vanilla Adam (https://arxiv.org/abs/1412.6980).
export Adam, step!

# Struct containing all necessary info
mutable struct Adam
  theta::AbstractArray{Float64} # Parameter array
  loss::Function                # Loss function
  grad::Function                # Gradient function
  m::AbstractArray{Float64}     # First moment
  v::AbstractArray{Float64}     # Second moment
  b1::Float64                   # Exp. decay first moment
  b2::Float64                   # Exp. decay second moment
  a::Float64                    # Step size
  eps::Float64                  # Epsilon for stability
  t::Int                        # Time step (iteration)
end

# Outer constructor
function Adam(theta::AbstractArray{Float64}, loss::Function, grad::Function)
  m   = zeros(size(theta))
  v   = zeros(size(theta))
  b1  = 0.9
  b2  = 0.999
  a   = 0.1
  eps = 1e-8
  t   = 0
  Adam(theta, loss, grad, m, v, b1, b2, a, eps, t)
end

# Step function with optional keyword arguments for the data passed to grad()
function step!(opt::Adam; data...)
  opt.t += 1
  gt    = opt.grad(opt.theta; data...)
  opt.m = opt.b1 .* opt.m + (1 - opt.b1) .* gt
  opt.v = opt.b2 .* opt.v + (1 - opt.b2) .* gt .^ 2
  mhat = opt.m ./ (1 - opt.b1^opt.t)
  vhat = opt.v ./ (1 - opt.b2^opt.t)
  opt.theta -= opt.a .* (mhat ./ (sqrt.(vhat) .+ opt.eps))
end

end
