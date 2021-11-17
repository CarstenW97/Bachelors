fig = Figure(fontsize = 28, backgroundcolor = :transparent);
# fig = Figure(font = assetpath("fonts", fontpath), fontsize=24);
ax = Axis(fig[1, 1]);


b = 0.05 # [0, 0.01, 0.05, 0.1, 0.15]
d_f = 2 # [2,10,20,30,40,50,75,100]

dft = filter([:driving_force, :b] => (x, y) -> x == d_f && y == b, df)

idxs = sortperm(dft.:S_ext)
xs = dft.:S_ext[idxs]

prot_b = dft.:b[idxs]
prot_tw = dft.:b[idxs] + dft.:tw[idxs]
prot_ts = dft.:b[idxs] + dft.:tw[idxs] + dft.:ts[idxs]
prot_l = dft.:b[idxs] + dft.:tw[idxs] + dft.:ts[idxs] + dft.:l[idxs]
prot_r = dft.:b[idxs] + dft.:tw[idxs] + dft.:ts[idxs] + dft.:l[idxs] + dft.:r[idxs]
prot_a = dft.:b[idxs] + dft.:tw[idxs] + dft.:ts[idxs] + dft.:l[idxs] + dft.:r + dft.:a[idxs]
prot_c =
    dft.:b[idxs] +
    dft.:tw[idxs] +
    dft.:ts[idxs] +
    dft.:l[idxs] +
    dft.:r +
    dft.:a[idxs] +
    dft.:c[idxs]

band!(
    ax,
    xs,
    zeros(length(prot_b)),
    prot_b ./ prot_c,
    color = colscheme[1],
    label = "Transporter",
)
band!(
    ax,
    xs,
    prot_b ./ prot_c,
    prot_tw ./ prot_c,
    color = colscheme[2],
    label = "Metabolic enzyme",
)
band!(
    ax,
    xs,
    prot_tw ./ prot_c,
    prot_ts ./ prot_c,
    color = colscheme[3],
    label = "Ribosomes",
)
band!(
    ax,
    xs,
    prot_ts ./ prot_c,
    prot_l ./ prot_c,
    color = colscheme[4],
    label = "Lipid synthesis enzyme",
)


#hidexdecorations!(ax, ticks = false, ticklabels = false, label=false)
#hideydecorations!(ax, ticks = false, ticklabels = false, label=false)
ax.xlabel = "External substrate concentration"
ax.ylabel = "Proteome fraction"
fig[1, 2] = Legend(fig, ax, "Protein", unique = true, framevisible = false)


fig
