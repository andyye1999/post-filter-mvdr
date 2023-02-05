./apply-delay-and-sum --tdoa-window=4000 --beam-window=4000 --margin=16 \
    test1.wav test1.delay_and_sum.wav

./apply-gsc --num-k=128 --alpha=0.01 \
    test1.wav test1.gsc.wav

./apply-mvdr --frame-len=0.025 --frame-shift=0.01 --fft-point=512 \
    --energy-thresh=1.5e-7 \
    --sil-to-speech-trigger=3 \
    --speech-to-sil-trigger=10 \
    test1.wav test1.mvdr.wav