import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  server: {
    port: 4000,
    proxy: {
      // In dev, forward API calls to the backend without CORS issues
      '/test-soft-powers': 'http://localhost:8200',
      '/analyze':          'http://localhost:8200',
      '/status':           'http://localhost:8200',
      '/results':          'http://localhost:8200',
      '/files':            'http://localhost:8200',
    },
  },
})
