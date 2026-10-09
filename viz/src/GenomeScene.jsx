import { useRef, useMemo, useEffect, useState } from 'react'
import { Canvas, useFrame, useThree } from '@react-three/fiber'
import { OrbitControls, Text, Line } from '@react-three/drei'
import * as THREE from 'three'

// Chromosome lengths in Mb (GRCh38 approximate)
const CHR_LENGTHS = [
  248,242,198,190,181,170,159,145,138,133,
  135,133,114,107,101,90,83,80,58,64,46,50
];
const TOTAL_MB = CHR_LENGTHS.reduce((a,b)=>a+b,0);

const PATHWAY_COLORS = {
  'RAS/MAPK': '#58a6ff',
  'RTK/Growth': '#3fb950',
  'Cell Cycle': '#f85149',
  'Chromatin': '#bc8cff',
  'DNA Repair': '#e3b341',
  'PI3K/AKT': '#ff7b72',
  'Wnt': '#39d0d8',
};

function ChromosomeRing({ index, length, variants, selectedGene, onVariantClick }) {
  const meshRef = useRef()
  const ringRadius = 7
  const angle = (index / 22) * Math.PI * 2
  const cx = Math.cos(angle) * ringRadius
  const cz = Math.sin(angle) * ringRadius
  const chrHeight = (length / 248) * 3.2

  useFrame((state) => {
    if (meshRef.current) {
      meshRef.current.rotation.y = state.clock.elapsedTime * 0.05
    }
  })

  const isActive = variants.length > 0

  return (
    <group position={[cx, 0, cz]}>
      {/* Chromosome body */}
      <mesh ref={meshRef}>
        <cylinderGeometry args={[0.12, 0.12, chrHeight, 8]} />
        <meshStandardMaterial
          color={isActive ? '#1a2744' : '#0d1117'}
          emissive={isActive ? '#0a1529' : '#000'}
          transparent opacity={0.9}
        />
      </mesh>
      {/* Chromosome glow if has variants */}
      {isActive && (
        <mesh>
          <cylinderGeometry args={[0.18, 0.18, chrHeight, 8]} />
          <meshStandardMaterial
            color="#58a6ff" transparent opacity={0.04}
            wireframe={false}
          />
        </mesh>
      )}
      {/* Chr label */}
      <Text
        position={[0, chrHeight / 2 + 0.3, 0]}
        fontSize={0.18}
        color={isActive ? '#7d9dcc' : '#2a3347'}
        anchorX="center" anchorY="bottom"
      >
        {index + 1}
      </Text>
      {/* Variant dots */}
      {variants.map((v, i) => {
        const fracAlong = (v.pos - (v.pos % 1000000)) / (CHR_LENGTHS[index] * 1e6)
        const yPos = (fracAlong - 0.5) * chrHeight
        const isSelected = selectedGene === v.gene
        const r = isSelected ? 0.22 : 0.14
        const color = PATHWAY_COLORS[v.pathway] || v.color || '#58a6ff'
        return (
          <group key={i} position={[0, yPos, 0]}>
            <mesh
              onClick={(e) => { e.stopPropagation(); onVariantClick(v) }}
              onPointerOver={(e) => { e.stopPropagation(); document.body.style.cursor='pointer' }}
              onPointerOut={() => { document.body.style.cursor='auto' }}
            >
              <sphereGeometry args={[r, 8, 8]} />
              <meshStandardMaterial
                color={color}
                emissive={color}
                emissiveIntensity={isSelected ? 2.5 : 1.2}
              />
            </mesh>
            {isSelected && (
              <mesh>
                <sphereGeometry args={[0.30, 8, 8]} />
                <meshStandardMaterial color={color} transparent opacity={0.15} />
              </mesh>
            )}
          </group>
        )
      })}
    </group>
  )
}

function CenterDNA({ processing }) {
  const groupRef = useRef()
  const points = useMemo(() => {
    const pts = []
    for (let i = 0; i <= 60; i++) {
      const t = i / 60
      const angle = t * Math.PI * 4
      pts.push(new THREE.Vector3(Math.cos(angle) * 0.6, t * 4 - 2, Math.sin(angle) * 0.6))
    }
    return pts
  }, [])
  const points2 = useMemo(() => points.map(p =>
    new THREE.Vector3(-p.x, p.y, -p.z)
  ), [points])

  useFrame((s) => {
    if (groupRef.current) {
      groupRef.current.rotation.y = s.clock.elapsedTime * (processing ? 1.5 : 0.3)
    }
  })

  const rungs = useMemo(() => {
    const r = []
    for (let i = 2; i < 60; i += 4) {
      r.push([points[i], points2[i]])
    }
    return r
  }, [points, points2])

  return (
    <group ref={groupRef} position={[0, 0, 0]}>
      <Line points={points} color={processing ? '#3fb950' : '#1a3a5c'} lineWidth={2} />
      <Line points={points2} color={processing ? '#58a6ff' : '#1a3a5c'} lineWidth={2} />
      {rungs.map((r, i) => (
        <Line key={i} points={r}
          color={processing ? '#bc8cff' : '#0d2137'} lineWidth={1} opacity={0.6} transparent />
      ))}
    </group>
  )
}

function ParticleField({ active }) {
  const ref = useRef()
  const count = 300
  const positions = useMemo(() => {
    const pos = new Float32Array(count * 3)
    for (let i = 0; i < count; i++) {
      pos[i*3]   = (Math.random()-0.5) * 24
      pos[i*3+1] = (Math.random()-0.5) * 12
      pos[i*3+2] = (Math.random()-0.5) * 24
    }
    return pos
  }, [])

  useFrame((s) => {
    if (ref.current) {
      ref.current.rotation.y = s.clock.elapsedTime * 0.02
    }
  })

  return (
    <points ref={ref}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" args={[positions, 3]} />
      </bufferGeometry>
      <pointsMaterial size={0.03} color={active ? '#2a4a7f' : '#1a2535'} transparent opacity={0.6} />
    </points>
  )
}

function Scene({ variants, selectedGene, onVariantClick, processing }) {
  // Group variants by chromosome
  const byChr = useMemo(() => {
    const map = {}
    variants.forEach(v => {
      const c = v.chr - 1
      if (!map[c]) map[c] = []
      map[c].push(v)
    })
    return map
  }, [variants])

  return (
    <>
      <color attach="background" args={['#060b14']} />
      <ambientLight intensity={0.5} />
      <pointLight position={[0, 10, 0]} intensity={3} color="#3060a0" />
      <pointLight position={[0, -5, 0]} intensity={1.5} color="#204060" />
      <pointLight position={[10, 0, 10]} intensity={1.5} color="#102030" />
      <ParticleField active={processing} />
      <CenterDNA processing={processing} />
      {CHR_LENGTHS.map((len, i) => (
        <ChromosomeRing
          key={i} index={i} length={len}
          variants={byChr[i] || []}
          selectedGene={selectedGene}
          onVariantClick={onVariantClick}
        />
      ))}
      <OrbitControls
        enablePan={true} enableZoom={true} enableRotate={true}
        minDistance={5} maxDistance={35}
        autoRotate={!processing} autoRotateSpeed={0.4}
      />
    </>
  )
}

export default function GenomeScene({ variants, selectedGene, onVariantClick, processing }) {
  return (
    <Canvas
      camera={{ position: [0, 8, 18], fov: 55 }}
      style={{ width: '100%', height: '100%' }}
      dpr={[1, 1.5]}
    >
      <Scene
        variants={variants}
        selectedGene={selectedGene}
        onVariantClick={onVariantClick}
        processing={processing}
      />
    </Canvas>
  )
}
