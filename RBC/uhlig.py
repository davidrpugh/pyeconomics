class RBC:
    """
    
    Class RBC is used to solve and simulate a RBC models a la methods detailed in Uhlig (1997)!
    
    """
    
    def __init__(self, coefs=None, P=None, Q=None, R=None, S=None, x=None, y=None, z=None, e=None, phi=None):
        """
        
        Attributes: 
            
            1) coefs is a list containing the arrays (matrices) A, ..., N collecting the coefficients 
               of the model (see section 6.2 of Uhlig (1997) for more details).
            2) If supplied (i.e., you have already solved the model by hand using method of 
               undetermined coefficients!), then P, Q, R, and S are matrices describing the 
               recursive law of motion.   
            3) x, is a vector of initial conditions for the endogenous state variables
            4) y, is a vector of initial conditions for the other endogenous variables of the model
            5) z, is a vector of initial conditions for the exogenous state variables (i.e., shocks)
            6) e, is a vector of initial conditions for the innovations driving exogenous shocks
            7) phi, is a function returning either a vector of specified shocks, or a vector of draws 
               from the distributions of innovations of the z variables.
        
        """
        self.coefs, self.P, self.Q, self.R, self.S, self.x, self.y, self.z, self.e, self.phi = coefs, P, Q, R, S, x, y, z, e, phi
    
    def pre_solve(self):
        """
        
        Method pre_solve() computes the matrices PSI, GAMMA, THETA, XI, and 
        DELTA from the supplied coefficient matrices A,...,N.  The matrices 
        A,...,N can be used to define matrix quadratic equation in the 
        recursive law of motion for the endogenous state variable(s), P:
            
            PSI * P**2 - GAMMA * P - THETA = 0
        
        In general, this quadratic equation is solved by, first, using the 
        matrices PSI, GAMMA, and THETA to construct the matrices XI and DELTA, 
        and, second, constructing P in terms of stable eigenvalues and eigen-
        vectors that are solutions to the following generalized eigen-value 
        problem:
            
            XI * x = lambda * DELTA * x
            
        Returns a list containing [k, l, m, n, q, PSI, GAMMA, THETA, XI, DELTA]
        
        """
        
        # Unpack the list of coefficient matrices
        A, B, C, D = self.coefs[0], self.coefs[1], self.coefs[2], self.coefs[3]
        F, G, H, J, K, L, M = self.coefs[4], self.coefs[5], self.coefs[6], self.coefs[7], self.coefs[8], self.coefs[9], self.coefs[10]
        N = self.coefs[11]
        
        # The number of columns (or rows) of the matrix N equals the number of shocks
        k = np.size(N, axis=0)
        
        # The number of columns in A equals the number of endogenous state variables
        m = np.size(A, axis=1)
        
        # The number of columns in C equals the number of endogenous non-state variables
        n = np.size(C, axis=1)
        
        # The number of rows in C is the number of non-expectational equations
        l = np.size(C, axis=0)
        
        # The number of rows of F equals the number of expectational equations
        q = np.size(F, axis=0)
                 
        # C_plus is the Moore-Penrose (or pseudo-inverse) of the matrix C 
        C_plus = linalg.pinv(C)
        
        def null(X, tol=1e-15):
            """
            
            Calculates a basis for the null space of X using SVD.
            
            Inputs:
                1) an m x n matrix, X
                2) tol, specifies how close to zero the singular values
                   should be in order to be considered equal to zero.
            
            Returns:
                1) matrix whose rows form a basis for the null space of X (i.e.,
                   those rows of V_trans whose singular values are zero!)
                
            """ 
            U, Sigma, V_trans = np.linalg.svd(X)
            null_mask = (Sigma <= tol)
            if np.sum(null_mask) == 0:
                return np.zeros(shape=(np.size(null_mask), 1))
            else:
                null_space = np.compress(null_mask, V_trans, axis=0)
                return null_space
        
        # C_0 is a matrix whose ROWS form a basis for the null space of C'
        C_0 = null(C.T)
        
        # Define the PSI, GAMMA, and THETA matrices...
        if l == 0:
            """
            
            You are using the 'Brute force" approach from section 6.1 of Uhlig (1997)!
            The number of expectational equations, q, should equal the number of endogenous
            variables, n!
            
            """
            self.PSI = F
            self.GAMMA = -G
            self.THETA = -H
        elif n == l:
            """
            
            You are using the 'Sensitive' approach from section 6.2 of Uhlig (1997)!
            Your model is simple (i.e., has the same number of equations as endogenous 
            (non-state) variables and you get to use simpler formulas for PSI, GAMMA, 
            and THETA!
            
            Note that when l = n, C is invertible and thus C_plus = sp.pinv(C) = sp.inv(C)!
            
            """
            self.PSI = F - J.dot(C_plus).dot(A)
            self.GAMMA = J.dot(C_plus).dot(B) - G + K.dot(C_plus).dot(A)
            self.THETA = K.dot(C_plus).dot(B) - H
        else:
            """
            
            You are using the 'Sensitive' approach from section 6.2 of Uhlig (1997)!
            However, you have more equations than (non-state) endogenous variables (i.e.,
            l > n).
            
            """
            self.PSI = np.vstack(tup=(np.zeros(shape=(l - n, m)), F - J.dot(C_plus).dot(A))) 
            self.GAMMA = np.vstack(tup=(C_0.dot(A), J.dot(C_plus).dot(B) - G + K.dot(C_plus).dot(A)))
            self.THETA = np.vstack(tup=(C_0.dot(B), K.dot(C_plus).dot(B) - H))
            
        # Construct the XI and DELTA matrices used in solving the matrix quadratic equation 
        self.XI = np.vstack(tup=(np.hstack(tup=(self.GAMMA, self.THETA)), np.hstack(tup=(np.eye(m), np.zeros(shape=(m, m))))))
        self.DELTA = np.vstack(tup=(np.hstack(tup=(self.PSI, np.zeros(shape=(m, m)))), np.hstack(tup=(np.zeros(shape=(m, m)), np.eye(m)))))
        
        return [k, l, m, n, q, self.PSI, self.GAMMA, self.THETA, self.XI, self.DELTA]
     
    def solve_P(self, method='eigen'):
        """
        
        Method solves for the recursive law of motion for the endogenous state 
        variable(s) of the model model (i.e., transition matrix P) using either:
            
            1) method = 'eigen': Generalized eigenvalue routine detailed in 
               Uhlig (1997).  Method returns a list containing: 
                   
                   [XI_eigenvalues, XI_eigenvectors, OMEGA, LAMBDA, P]
                   
               First two elements of this list are self explanatory.  OMEGA and LAMBDA
               are the m eigenvalues and corresponding eigenvalues selected for use in 
               calculating P. 
            2) method = 'QZ': QZ decomposition due to Sims (1996)
        
        """
        k, l, m, n, q, PSI, GAMMA, THETA, XI, DELTA = self.pre_solve()
        
        if method == 'eigen':
            """
            
            Solve for recursive law of motion via generalized eigenvalue decomposition as
            detailed in Uhlig (1997). In order to have a stable solution, P, we need to find 
            m out of the 2 * m possible generalized eigenvalues and their associated 
            eigenvectors such that |lambda_i| < 1 for all i = 1,...,m.
            
            Ideally, there will be exactly m such eigenvalues and eigenvectors (i.e., your
            model is 'saddle-point stable'). However, this need not always be the case.  
            Particular, if your model has multiple equilibria and/or multiple steady-states!
            
            """
            # Solve the generalized eigenvalue problem
            self.XI_eigenvalues, self.XI_eigenvectors = linalg.eig(XI, DELTA)
            # Ignore imaginary part of eigenvalues/vectors if imaginary part is "close" to zero
            self.XI_eigenvalues = np.real_if_close(self.XI_eigenvalues, tol=100)
            self.XI_eigenvectors =  np.real_if_close(self.XI_eigenvectors, tol=100)
             
            if np.linalg.matrix_rank(self.XI_eigenvectors) < m:
                # Need at least m linearly independent eigenvectors to solve for P! 
                return 'Sorry! XI is not diagonalizable! Try using method = QZ instead!'  
            else:
                # Identify stable eigenvalues and their associated eigenvectors
                XI_mask = (np.abs(self.XI_eigenvalues) < 1)
                if np.sum(XI_mask) > m:
                    return 'Too many stable eigenvalues!'
                self.LAMBDA = np.diag(np.compress(XI_mask, self.XI_eigenvalues, axis=0))
                self.OMEGA = np.compress(XI_mask, np.compress(XI_mask, self.XI_eigenvectors, axis=0), axis = 1)
                if np.linalg.matrix_rank(self.OMEGA) < m:
                    return 'Sorry! OMEGA is not invertible! Try using method = QZ instead!'
                # Finally! Calculate P
                self.P = self.OMEGA.dot(self.LAMBDA).dot(linalg.inv(self.OMEGA))
                
            #return [self.XI_eigenvalues, self.XI_eigenvectors, self.OMEGA, self.LAMBDA]
        
        if method == 'QZ':
            """
            
            Implementation of the QZ decomposition a la Sims (1996)
            
            """
            # TO BE WRITTEN!
        
        
    def solve_QRS(self):
        """
        
        Calculates the matrices Q, R, and S which (combined with P) fully 
        describe the recursive law of motion for the system.
        
        """
        
        # Unpack the list of coefficient matrices
        A, B, C, D = self.coefs[0], self.coefs[1], self.coefs[2], self.coefs[3]
        F, G, H, J, K, L, M = self.coefs[4], self.coefs[5], self.coefs[6], self.coefs[7], self.coefs[8], self.coefs[9], self.coefs[10]
        N = self.coefs[11]
        
        # The number of columns (or rows) of the matrix N equals the number of shocks
        k = np.size(N, axis=0)
        
        # The number of columns in A equals the number of endogenous state variables
        m = np.size(A, axis=1)
        
        # The number of columns in C equals the number of endogenous non-state variables
        n = np.size(C, axis=1)
        
        # The number of rows in C is the number of non-expectational equations
        l = np.size(C, axis=0)
        
        # The number of rows of F equals the number of expectational equations
        q = np.size(F, axis=0)
                 
        # C_plus is the Moore-Penrose (or pseudo-inverse) of the matrix C 
        C_plus = linalg.pinv(C)
        
        def vec(X):
            """
            
            Function returns the vectorization of the array X by first 
            flattening the array to remove the original dimensions, and 
            then creating a single column vector out of the flattened
            array.
            
            """
            rows_X, cols_X = np.shape(X)[0], np.shape(X)[1]
            tmp_X = X.flatten()
            vec_X = np.reshape(tmp_X, newshape=(rows_X * cols_X, 1))
            return vec_X
            
        if l == 0:
            """
            
            You are using the 'Brute force" approach from section 6.1 of Uhlig (1997)!
            The number of expectational equations, q, should equal the number of endogenous
            variables, n!
            
            Using the brute force method means that the R and S matrices are not defined.
            
            """
            # Q is a bit complicated. First, the LHS of 6.5...
            V = linalg.kron(N.T, F) + linalg.kron(np.eye(k), F.dot(self.P) + G)
            
            #...and the the RHS of 6.5
            vec_U = -vec(L.dot(N) + M)
            
            if np.linalg.matrix_rank(V) < k * m:
                print 'V is not invertible! Cannot solve for Q!'
            else:
                vec_Q = linalg.inv(V).dot(vec_U)

        elif n == l:
            """
            
            You are using the 'Sensitive' approach from section 6.2 of Uhlig (1997)!
            Your model is simple (i.e., has the same number of equations as endogenous 
            (non-state) variables and you get to use simpler formulas for Q, R, and S!
            
            Note that when l = n, C is invertible and thus C_plus = sp.pinv(C) = sp.inv(C)!
            
            """
            # Given P, R is simple!
            self.R = -C_plus.dot(A.dot(self.P) + B)
            
            # Q is more complicated. First, the LHS of 6.19...
            V_11 = linalg.kron(N.T, F - J.dot(C_plus).dot(A))
            V_12 = linalg.kron(np.eye(k), J.dot(self.R) + F.dot(self.P) + G - K.dot(C_plus).dot(A))
            V = V_11 + V_12
            
            #...and the the RHS of 6.19
            U = (J.dot(C_plus).dot(D) - L).dot(N) + K.dot(C_plus).dot(D) - M
            vec_U = vec(U)
            
            if np.linalg.matrix_rank(V) < k * m:
                print 'V is not invertible! Cannot solve for Q, R, and S!'
            else:
                vec_Q = linalg.inv(V).dot(vec_U)
                self.Q = np.reshape(vec_Q, newshape=(m, k))
  
            # Given Q, S is simple!
            self.S = -C_plus.dot(A.dot(self.Q) + D)

        else:
            """
            
            You are using the 'Sensitive' approach from section 6.2 of Uhlig (1997)!
            
            """
            # Given P, R is simple!
            self.R = -C_plus.dot(A.dot(self.P) + B)
            
            # Q is a bit more complicated. First, the LHS of 6.14...
            V_11 = linalg.kron(np.eye(k), A)
            V_12 = linalg.kron(np.eye(k), C)
            V_21 = linalg.kron(N.T, F) + linalg.kron(np.eye(k), F) + linalg.kron(np.eye(k), J.dot(self.R) + F.dot(self.P) + G)
            V_22 = linalg.kron(N.T, J) + linalg.kron(np.eye(k), K)
            V_1 = np.hstack(tup=(V_11, V_12))
            V_2 = np.hstack(tup=(V_21, V_22))
            V = np.vstack(tup=(V_1, V_2))
            
            #...and the the RHS of 6.14
            vec_U = -np.vstack(tup=(vec(D), vec(L.dot(N) + M)))
            
            if np.linalg.matrix_rank(V) < k * (m + n):
                print 'V is not invertible! Cannot solve for Q, R, and S!'
            else:
                vec_Q = linalg.inv(V).dot(vec_U)[0:(k * m -1), 0]
                vec_S = linalg.inv(V).dot(vec_U)[(k * m):(k * (m + n) - 1), 0]
                self.Q = np.reshape(vec_Q, newshape=(m, k))
                self.S = np.reshape(vec_Q, newshape=(n, k))
        
    def update(self, mode, shock):
        """
        
        Uses the matrices describing the recursive law of motion (i.e., 
        P, Q, R, and S) to updates the endogenous variables according to: 
            
            x_t = P * x_{t-1} + Q * z_t
            y_t = R * x_{t-1} + S * z_t
        
        """
        N = self.coefs[11]
        self.z = N.dot(self.z) + self.phi(mode, shock)
        self.x = self.P.dot(self.x) + self.Q.dot(self.z)
        self.y = self.R.dot(self.x) + self.S.dot(self.z)
        
    def impulse_response(self, n, t, shock):
        """
        
        Calculates the impulse responses of the endogenous variables of the model
        to an exogenous shock.  
        
        Inputs:
            
            1) n: length of impulse response to calculate
            2) t: period in which the shock hits
            3) shock: vector of length k defining the size of the shocks
            
        Returns:
            
            1) a list containing the impulse response functions
               for the variables of the model 
        
        """
        N = self.coefs[11]
        
        # Create arrays to hold sample paths...
        path_x = np.zeros(shape=(np.shape(self.P)[0], n))
        path_y = np.zeros(shape=(np.shape(self.R)[0], n))
        path_z = np.zeros(shape=(np.shape(N)[0], n))
        path_e = np.zeros(shape=(np.shape(N)[0], n))
                
        # Shock is set to zero for n != t!
        for i in range(1, n):
            if i == t:
                self.update('IR', shock)
                path_x[:, t] = self.x.flatten()
                path_y[:, t] = self.y.flatten()
                path_z[:, t] = self.z.flatten()
                path_e[:, t] = self.e.flatten()         
            else:
                self.update('IR', np.zeros(shape=(2,1)))
                path_x[:, i] = self.x.flatten()
                path_y[:, i] = self.y.flatten()
                path_z[:, i] = self.z.flatten()
                path_e[:, i] = self.e.flatten()             
        
        return [path_x.T, path_y.T, path_z.T, path_e.T]

    def sample_path(self, n, seed):
        """
        
        Generate path of length n for both the endogenous variable and the 
        stochastic shock process.  
        
        Returns a list containing the simulated time paths.  Note that time
        paths are returned as column vectors!
        
        """
        # set the seed
        np.random.seed(int(seed))
        
        N = self.coefs[11]
        
        path_x = np.zeros(shape=(np.shape(self.P)[0], n))
        path_y = np.zeros(shape=(np.shape(self.R)[0], n))
        path_z = np.zeros(shape=(np.shape(N)[0], n))
        path_e = np.zeros(shape=(np.shape(N)[0], n))
     
        for t in range(n):
            path_x[:, t] = self.x.flatten()
            path_y[:, t] = self.y.flatten()
            path_z[:, t] = self.z.flatten()
            path_e[:, t] = self.e.flatten()
            self.update(mode='simul', shock=None)
        
        return [path_x.T, path_y.T, path_z.T, path_e.T]
