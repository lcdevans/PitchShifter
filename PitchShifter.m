clc
clf
clear all

% Read an audio file
soundFile = 'F:\Fall 2020\ECE 484\sample1.wav';
[input1, Fs] = audioread(soundFile);
%sound(input, Fs);

numChannels = size(input1,2);

% Set up parameters
hopSizeA = 32;
frameSize = 2048;
figure(1)
plot(abs(fft(input1, Fs)))
title("Frequency Magnitude Plot for 440 Hz Input");
xlabel("Frequency (Hz - 1)");
ylabel("Magnitude")

Ratio = 1.06;
[num, dem] = rat(Ratio);
hopSizeS = int64(Ratio * hopSizeA);

numHops = ceil(size(input1,1)/hopSizeA);
magsA = zeros(numHops, frameSize, numChannels);
phasesA = zeros(numHops, frameSize, numChannels);

magsS = zeros(numHops, frameSize, numChannels);
phasesS = zeros(numHops, frameSize, numChannels);
framesS = zeros(numHops, frameSize, numChannels);

% Pad the input array
padArray = zeros(frameSize, numChannels);
input = [input1;padArray];
output = zeros(3*size(input, 1), numChannels);

% Hann window
w = 0.5*(1-cos(2*pi*(0:n-1)'/(n)));
wS = w;

% Stacking windows with hop ratio < 0.25 will become large and start
% clipping.
% Pre-calculate this amount so we can scale the window accordingly.
testWindowA = zeros(2 * frameSize, 1);
numHopPerFrameA = frameSize/hopSizeA;
for hop = 1:numHopPerFrameA
    testStartA = (hop - 1) * hopSizeA + 1;
    testEndA = testStartA + frameSize - 1;
    testWindowA(testStartA:testEndA) = testWindowA(testStartA:testEndA) + w;
end
maxA = max(testWindowA);
w = w/maxA;

testWindowS = zeros(2 * frameSize, 1);
numHopPerFrameS = frameSize/hopSizeS;
for hop = 1:numHopPerFrameS
    testStartS = (hop - 1) * hopSizeS + 1;
    testEndS = testStartS + frameSize - 1;
    testWindowS(testStartS:testEndS) = testWindowS(testStartS:testEndS) + wS;
end
maxS = max(testWindowS);
wS = wS/maxS;


% Perform STFT Analysis step
for channel = 1:numChannels
    for hopA = 1:numHops
        frameStartA = (hopA - 1) * hopSizeA + 1;
        frameEndA = frameStartA + frameSize - 1;
        window = w .* input(frameStartA:frameEndA, channel);
        window = circshift(window, frameSize/2);
        framefft = fft(window);
        magsA(hopA, :, channel) = abs(framefft).';
        phasesA(hopA, :, channel) = atan2(imag(framefft), real(framefft)).';
    end
end

omegakh = 2 * pi * (0:frameSize-1) * hopSizeA / frameSize;
%Do Phase stuff
for channel = 1:numChannels
    for hop = 1:numHops
        
      if hop == 1
          magsS(hop, :, channel) = magsA(hop, :, channel);
          phasesS(hop, :, channel) = phasesA(hop, :, channel);
      else
          magsS(hop, :, channel) = magsA(hop, :, channel);
          phaseDiff = phasesA(hop, :, channel) - phasesA(hop - 1, :, channel);
          test1 = phaseDiff - omegakh;
          test1 = test1 - round(test1/(2*pi))*2*pi;
          test1 = (test1 + omegakh) * Ratio;
          phasesS(hop, :, channel) = test1 + phasesS(hop-1, :, channel);
      end
      
      framesS(hop, :, channel) = magsS(hop, :, channel) .* cos(phasesS(hop, :, channel)) + 1i * magsS(hop, :, channel) .* sin(phasesS(hop, :, channel));
    end
end

% Synthesis step
for channel = 1:numChannels
    for hop = 1:numHops
        frameStartS = (hop - 1) * hopSizeS + 1;
        frameEndS = frameStartS + frameSize - 1;
        frameSifft = real(ifft(framesS(hop, :, channel)));
        windowS = circshift(frameSifft, frameSize/2);
                
        output(frameStartS:frameEndS, channel) = output(frameStartS:frameEndS, channel) + windowS.';
    end
end
output1 = resample(output, dem, num);
figure(2)
plot(abs(fft(output1, Fs)))
title("Frequency Magnitude Plot for Hop Size = 32, Frame Size = 2048, Semitone Up");
xlabel("Frequency (Hz - 1)");
ylabel("Magnitude")
sound(output1, Fs)

%memes = input1;
%memes(:,1) = output1(1:size(input1,1), 2);
%sound(memes, Fs)
