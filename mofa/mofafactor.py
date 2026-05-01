import numpy as np

class MOFAStyleModel:
    def __init__(
        self,
        k,
        params,
        iters=100,
        tol=1e-4,
        reg=None,
        random_state=0,
        normalize_factors=True,
        verbose=True,
        print_iter=10
    ):
        self.k = k
        self.params = params
        self.iters = iters
        self.tol = tol
        self.reg = reg
        self.random_state = random_state
        self.normalize_factors = normalize_factors
        self.verbose = verbose
        self.print_iter = print_iter

        self.Z = None          # (n x k)
        self.W = {}            # dict of (d_v x k)

    # -------------------------
    # Loss
    # -------------------------
    def loss(self, Xs):
        loss = 0
        for v, X in Xs.items():
            Wv = self.W[v]
            lam = self.params[f"lamda_{v}"]
            loss += lam * np.linalg.norm(X - self.Z @ Wv.T)**2 / 2

        # regularization on Z
        mu_1 = self.params.get("mu_1", 0)
        mu_2 = self.params.get("mu_2", 0)

        loss += mu_1 * np.sum(np.abs(self.Z))
        loss += mu_2 * np.linalg.norm(self.Z)**2 / 2

        if self.reg is not None:
            D = self.reg["D_lap"]
            A = self.reg["A_lap"]
            mu_lap = self.params.get("mu_lap", 0)
            loss += mu_lap * np.trace(self.Z.T @ (D - A) @ self.Z) / 2

        return loss

    # -------------------------
    # Fit
    # -------------------------
    def fit(self, Xs):
        np.random.seed(self.random_state)

        self.views = list(Xs.keys())
        n = list(Xs.values())[0].shape[0]

        # init Z
        self.Z = np.random.rand(n, self.k)

        # init W
        for v, X in Xs.items():
            d = X.shape[1]
            self.W[v] = np.random.rand(d, self.k)

        trace = []
        obj = 1e10

        for i in range(self.iters):

            # -------------------------
            # Update Z (shared factors)
            # -------------------------
            num = 0
            denom = 0

            for v, X in Xs.items():
                Wv = self.W[v]
                lam = self.params[f"lamda_{v}"]

                num += lam * (X @ Wv)
                denom += lam * (self.Z @ (Wv.T @ Wv))

            # regularization
            mu_1 = self.params.get("mu_1", 0)
            mu_2 = self.params.get("mu_2", 0)

            denom += mu_1 + mu_2 * self.Z

            if self.reg is not None:
                A = self.reg["A_lap"]
                D = self.reg["D_lap"]
                mu_lap = self.params.get("mu_lap", 0)

                num += mu_lap * (A @ self.Z)
                denom += mu_lap * (D @ self.Z)

            self.Z *= num / (denom + 1e-10)

            if self.normalize_factors:
                self.Z /= self.Z.sum(axis=0, keepdims=True)

            # -------------------------
            # Update W (per view)
            # -------------------------
            for v, X in Xs.items():
                Wv = self.W[v]

                self.W[v] *= (X.T @ self.Z) / (Wv @ (self.Z.T @ self.Z) + 1e-10)

            # -------------------------
            # Check convergence
            # -------------------------
            if i % self.print_iter == 0:
                obj_new = self.loss(Xs)
                trace.append(obj_new)

                if self.verbose:
                    print(f"Iter {i}, loss={obj_new:.4f}")

                if abs(obj - obj_new) / obj < self.tol:
                    break

                obj = obj_new

        return trace

    # -------------------------
    # Reconstruction
    # -------------------------
    def reconstruction(self):
        return {v: self.Z @ self.W[v].T for v in self.views}

    # -------------------------
    # Getters
    # -------------------------
    def get_factors(self):
        return self.Z

    def get_weights(self, view):
        return self.W[view]
