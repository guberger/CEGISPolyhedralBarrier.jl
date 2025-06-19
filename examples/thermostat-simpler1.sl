(set-logic LIA)

(synth-inv InvF ((x Int) (t Int) (p Bool) (q Bool)))

(define-fun PreF ((x Int) (t Int) (p Bool) (q Bool)) Bool
    (and (>= x 16) (<= x 24) (= t 0) (not q))
)

(define-fun TransF ((x Int) (t Int) (p Bool) (q Bool)
                    (x! Int) (t! Int) (p! Bool) (q! Bool)) Bool
    (and
        (or
            (and p (>= x 19) (= x! (- x 1)) p!)
            (and p (<= x 18) (= x! (+ x 1)) (not p!))
            (and (not p) (<= x 21) (= x! (+ x 1)) (not p!))
            (and (not p) (>= x 22) (= x! (- x 1)) p!)
        )
        (or
            (and (not q) (<= t 49) (= t! (+ t 1)) (not q!))
            (and (not q) (>= t 49) (= t! 50) q!)
            (and q (>= t 49) (= t! 50) q!)
        )
    )
)

(define-fun PostF ((x Int) (t Int) (p Bool) (q Bool)) Bool
    (and
        (> x 0) (< x 30)
        (or (not q) (and (>= x 17) (<= x 23)))
    )
)

(inv-constraint InvF PreF TransF PostF)

(check-synth)