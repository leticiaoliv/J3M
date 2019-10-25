function declupivot(A::Matrix, T::type; diagtol = 1e-12)
    n = size(A, 2)
    p = collect(1:n)
    L = diagm(0 -> ones(T, n)
    U = convert(Matrix(T), A)
    for j = 1:n-1
        pivo, k = abs(A[j,j]), j
        for i = j+1:n
            if abs(A[i,j]) > pivo
                pivo, k = abs(A[i,j]), i
            end
        end
        if pivo <= diagtol
            error("Matriz singular ou muito próxima de ser singular")
        end

        if k != j
            p[k], p[j] = p[j], p[k]
            A[[k;j],:] = A[[j;k],:]
        end
        ajj = A[j,j]
        for i = j+1:n
            mij = A[i,j] / ajj
            A[i,j] = mij
            A[i,j+1:n] -= mij * A[j,j+1:n]
        end
    end
    L = tril(A, -1) + I
    U = triu(A)
    P = I * p
    return P, L, U
end

function solvelu(A, b, T)
    P, L, U = declupivot(copy(A),T)
    y = zeros(T, n)
    x = zeros(T, n)
    y = L \ P * b
    x = U \ y
    return x
end

function declurefine(A, b)
    x = solvelu(A, b, BigFloat)
    r = b - A * x
    ϵ = sqrt(eps(BigFloat))
    while abs(r) > ϵ
        Δx = solvelu(A, r, BigFloat)
        x += Δx
        r = b - A * x
    end
    return x
end

A = BigFloat.rand(100,100)
b = BigFloat.rand(100)
T = Float16
