// extension to real hermitian (symmetric) matrices 
numeric.eigh = function(A, maxiter) {
    if ('y' in A) { return numeric.jacobi_hermitian(A, maxiter); }
    return numeric.jacobi_real_symmetric(A, maxiter)
}


numeric.jacobinorm = function(A, B) {
    // used in numeric.jacobi
    var n = A.length;
    var s = 0;

    if (!B) {
        for (var i = 0; i < n; i ++) 
        {
            for (var j = i + 1; j < n; j ++) 
            {
                s = s + Math.pow(A[i][j], 2);
            }
        }
    } else {
        for (var i = 0; i < n; i ++) 
        {
            for (var j = i + 1; j < n; j ++) 
            {
                s = s + Math.pow(A[i][j], 2) + Math.pow(B[i][j], 2);
            }
        }
    }

    return Math.sqrt(s);
};

numeric.jacobi_real_symmetric = function(Ain, maxiter, options) {
    // jacobi method with rutishauser improvements from 
    // Rutishauser, H. (1966). The Jacobi method for real symmetric matrices. 
    // Numerische Mathematik, 9(1), 1–10. doi:10.1007/BF02165223

    // returns object containing
    // e: {x : v} eigenvalues.
    // lambda : {x: d} eigenstates.
    // niter : number of iterations.
    // iterations : list of convergence factors for each step of the iteration.
    // nrot : number of rotations performed.

    if (!options) { options = {}; }

    var size = [Ain.length, Ain[0].length];
    if (size[0] != size[1])
    {
        throw 'jacobi : matrix must be square';
    }
    // remember use only symmetric real matrices.
    var n = size[0];
    
    var v = numeric.identity(n);

    if (!('clone' in options ) || options['clone']===true) { 
        var A = numeric.clone(Ain);
    } else {
        var A = Ain;
    }

    var iters = numeric.rep([maxiter], 0);
    var d = numeric.getDiag(A);
    var bw = numeric.clone(d);

    // zeros    
    var zw = numeric.rep([n], 0);
    
    // iteration parameters
    var iter = -1;
    var niter = 0;
    var nrot = 0;
    var tresh = 0;
    
    //  prealloc
    var h, g, gapq, term, termp, termq, theta, t, c, s, tau;
    while (iter < maxiter)
    {
        iter++;
        iters[iter] = numeric.jacobinorm(A)
        niter = iter;
        tresh = iters[iter]/(4 * n);
        
        if (tresh==0) {
            return {E: {x: v, y: numeric.rep([n, n], 0)}, lambda: {x: d}, iterations: iters, niter: niter, nrot: nrot}; 
        }
        
        for (var p = 0; p < n; p++)
        {
            for (var q = p + 1; q < n; q++)
            {
                gapq = 10*Math.abs(A[p][q]);
                termp = gapq + Math.abs(d[p]);
                termq = gapq + Math.abs(d[q]);
                if (iter > 4 && termp == Math.abs(d[p]) && termq == Math.abs(d[q]))
                {
                    // remove small elmts
                    A[p][q] = 0;
                }
                else
                {
                    if (iter > 4 ||  Math.abs(A[p][q]) >= tresh)
                    {
                        // console.log('p: ' + p + ', q:' + q)

                        // apply rotation
                        h = d[q] - d[p];
                        term = Math.abs(h) + gapq;
                        if (term == Math.abs(h))
                        {
                            t = A[p][q]/h;
                        }
                        else
                        {
                            theta = 0.5 * h / A[p][q];
                            t = 1/(Math.abs(theta) + Math.sqrt(1 + theta*theta));
                            if (theta < 0)
                            {
                                t = -t;
                            }
                        }

                        // console.log(numeric.clone(A));

                        c = 1/Math.sqrt(1 + t * t);
                        s = t * c;
                        tau = s/(1 + c);
                        h = t * A[p][q];
                        zw[p] = zw[p] - h;
                        zw[q] = zw[q] + h;
                        d[p] = d[p] - h;
                        d[q] = d[q] + h;
                        A[p][q] = 0;
                        // console.log('h: ' + h);
                        // console.log('t: ' + t);
                        // console.log('c: ' + c);
                        // console.log('s: ' + s);
                        // console.log('tau:' + tau);
                        // console.log('s*tau: ' + (s * tau));

                        // rotate and use upper tria only
                        for (var j = 0; j < p; j++)
                        {
                            g = A[j][p];
                            h = A[j][q];
                            A[j][p] = g - s * (h + g * tau);
                            A[j][q] = h + s * (g - h * tau);
                        }
                        for (var j = p + 1; j < q; j++)
                        {
                            g = A[p][j];
                            h = A[j][q];
                            A[p][j] = g - s * (h + g * tau);
                            A[j][q] = h + s * (g - h * tau);
                        } 
                        for (var j = q + 1; j < n; j++)
                        {
                            g = A[p][j];
                            h = A[q][j];
                            A[p][j] = g - s * (h + g * tau);
                            A[q][j] = h + s * (g - h * tau);
                        }
                        // eigenstates
                        for (var j = 0; j < n; j++)
                        {
                            g = v[p][j];
                            h = v[q][j];
                            v[p][j] = g - s * (h + g * tau);
                            v[q][j] = h + s * (g - h * tau);
                        }
                        nrot++;
                    }
                }
            }
        }
        bw = numeric.add(bw, zw);
        d = numeric.clone(bw);
        zw = numeric.rep([n], 0);
    }

    return {E: numeric.t(v, numeric.rep([n, n], 0)), lambda: numeric.t(d, numeric.rep([n], 0)), iterations: iters, niter: niter, nrot: nrot};
};

numeric.jacobi_hermitian = function(Ain, maxiter, options) {
    // jacobi method with rutishauser improvements from 
    // Rutishauser, H. (1966). The Jacobi method for real symmetric matrices. 
    // Numerische Mathematik, 9(1), 1–10. doi:10.1007/BF02165223

    // returns object containing
    // E: {x : v, y: v} eigenvalues.
    // lambda : {x: d} eigenstates.
    // niter : number of iterations.
    // iterations : list of convergence factors for each step of the iteration.
    // nrot : number of rotations performed.
    if (!options) { options = {}; }

    var size = [Ain.x.length, Ain.x[0].length];
    if (size[0] != size[1])
    {
        throw 'jacobi : matrix must be square';
    }
    // remember use only symmetric real matrices.
    var n = size[0];
    
    if (!('clone' in options ) || options['clone']===true) { 
        var Ax = numeric.clone(Ain.x);
        var Ay = numeric.clone(Ain.y);
    } else {
        var Ax = Ain.x;
        var Ay = Ain.y;
    }
    
    var Ux = numeric.identity(n);
    var Uy = numeric.rep([n, n], 0);
    
    //
    var iters = numeric.rep([maxiter], 0);
    var d = numeric.getDiag(Ax);
    var bw = numeric.clone(d);

    // zeros    
    var zw = numeric.rep([n], 0);
    
    // iteration parameters
    var iter = -1;
    var niter = 0;
    var nrot = 0;
    var tresh = 0;
    
    //  prealloc
    var Apq, h, g, gapq, term, termp, termq, theta, t, c, s, tau;

    function sq(x, y) { return x * x + y * y; }
    function reprod(A, B) { return A.x * B.x - A.y * B.y; }
    function reprod1(A, B) { return A.x * B.x + A.y * B.y; }
    function reprod2(A, B) { return A.x * B.x + A.y * B.y; }
    function improd(A, B) { return A.x * B.y + A.y * B.x; }
    function improd1(A, B) { return A.x * B.y - A.y * B.x; }
    function improd2(A, B) { return -A.x * B.y + A.y * B.x; }

    while (iter < maxiter)
    {
        iter++;
        iters[iter] = numeric.jacobinorm(Ax, Ay);
        niter = iter;
        tresh = iters[iter]/(4 * n);

        if (tresh==0) { // note this works because of eps definition
            return {E: {x: Ux, y: Uy}, lambda: {x: d}, iterations: iters, niter: niter, nrot: nrot}; 
        }
        
        for (var p = 0; p < n; p++)
        {
            for (var q = p + 1; q < n; q++)
            {
                Apq = {x: Ax[p][q], y: Ay[p][q], sq: sq(Ax[p][q], Ay[p][q])};

                gapq = 10 * Apq.sq;
                termp = gapq + Math.abs(d[p]);
                termq = gapq + Math.abs(d[q]);
                
                if (iter > 4 && termp == Math.abs(d[p]) && termq == Math.abs(d[q]))
                {
                    // remove small elmts
                    Ax[p][q] = 0; Ay[p][q] = 0;
                }
                else 
                {
                    if (iter > 4 || Apq.sq >= tresh)
                    {
                        // console.log('p: ' + p + ', q:' + q)
                        // apply rotation
                        h = .5 * (d[q] - d[p]);
                        term = Math.abs(h) + gapq;
                        if (term == Math.abs(h))
                        {
                            t = .5 / h;
                        }
                        else
                        {
                            t = 1/(Math.abs(h) + Math.sqrt(Apq.sq + h * h));
                            if (h < 0)
                            {
                                t = -t;
                            }
                        }

                        h = t * Apq.sq;

                        c = 1 / Math.sqrt(1 + h * t);
                        s = - t * c;
                        tau = s * Apq.sq / (c + 1)

                        // zw
                        zw[p] = zw[p] - h;
                        zw[q] = zw[q] + h;

                        // shift ev
                        d[p] = d[p] - h;
                        d[q] = d[q] + h;

                        Ax[p][q] = 0; Ay[p][q] = 0;

                        // rotate and use upper tria only
                        for (var j = 0; j < p; j++)
                        {
                            g = {x: Ax[j][p], y: Ay[j][p]};
                            h = {x: Ax[j][q], y: Ay[j][q]};

                            Ax[j][p] = g.x + s * (reprod1(Apq, h) - tau * g.x);
                            Ay[j][p] = g.y + s * (improd1(Apq, h) - tau * g.y);
                            Ax[j][q] = h.x - s * (reprod(Apq, g) + tau * h.x);
                            Ay[j][q] = h.y - s * (improd(Apq, g) + tau * h.y);
                        }
                        for (var j = p + 1; j < q; j++)
                        {
                            g = {x: Ax[p][j], y: Ay[p][j]};
                            h = {x: Ax[j][q], y: Ay[j][q]};

                            Ax[p][j] = g.x + s * (reprod1(h, Apq) - tau * g.x);
                            Ay[p][j] = g.y + s * (improd1(h, Apq) - tau * g.y);
                            Ax[j][q] = h.x - s * (reprod1(g, Apq) + tau * h.x);
                            Ay[j][q] = h.y - s * (improd1(g, Apq) + tau * h.y);
                        } 
                        for (var j = q + 1; j < n; j++)
                        {
                            g = {x: Ax[p][j], y: Ay[p][j]};
                            h = {x: Ax[q][j], y: Ay[q][j]};
                            
                            Ax[p][j] = g.x + s * (reprod(Apq, h) - tau * g.x);
                            Ay[p][j] = g.y + s * (improd(Apq, h) - tau * g.y);
                            Ax[q][j] = h.x - s * (reprod1(Apq, g) + tau * h.x);
                            Ay[q][j] = h.y - s * (improd1(Apq, g) + tau * h.y);
                        }
                        // eigenstates
                        for( var j = 0; j < n; j++ ) {
                            g = {x: Ux[p][j], y: Uy[p][j]};
                            h = {x: Ux[q][j], y: Uy[q][j]};
                          
                            // Ux[p][j] = x.x + s * (reprod1(Apq, y) - t * x.x);
                            // Uy[p][j] = x.y + s * (improd1(Apq, y) - t * x.y);
                            // Ux[q][j] = y.x - s * (reprod(Apq, x) + t * y.x);
                            // Uy[q][j] = y.y - s * (improd(Apq, x) + t * y.y);
                            
                            Ux[p][j] = g.x + s * (reprod(Apq, h) - tau * g.x);
                            Uy[p][j] = g.y + s * (improd(Apq, h) - tau * g.y);
                            Ux[q][j] = h.x - s * (reprod1(Apq, g) + tau * h.x);
                            Uy[q][j] = h.y - s * (improd1(Apq, g) + tau * h.y);
                        }
                        nrot += 1;
                    }
                }
            }

        }
        bw = numeric.add(bw, zw);
        d = numeric.clone(bw);
        zw = numeric.rep([n], 0);
    }

    return {E: numeric.t(Ux, Uy), lambda: numeric.t(d, numeric.rep([n], 0)), iterations: iters, niter: niter, nrot: nrot}; 
};